#include <iostream>
#include <atomic>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <unordered_map>
#include <utility>
#include <optional>
#include <vector>
#include <stdexcept>
#include <filesystem>

#include <cinttypes>

#include "mpreal/mpreal.h"

#include "const_e.h"


#include "phys/units/io.hpp"
#include "phys/units/quantity.hpp"
#include "phys/units/other_units.hpp"


#define VERSION "0.3.0"
#define YEAR "2024"


using mpfr::mpreal;
using phys::units::quantity;

std::int32_t digits_prec = 0;
std::int32_t max_int_constants = 0;
std::int32_t max_expr_size = 1;

#define ERR_EXIT(A, ...) { \
    std::fprintf(stderr, "error: file " __FILE__ ":%i in %s(): ", __LINE__, __func__); \
    std::fprintf(stderr, __VA_ARGS__); \
    std::fputc('\n', stderr); \
    std::exit(static_cast<int>(A)); \
}

static constexpr const char CONSTANTS_FILENAME[] = "constants.conf";

static constexpr const char SAVE_AST_DIR[] = "save";

enum struct err_t : std::uint32_t {
    dimension_add, dimension_dim_exp, dimension_nonint_exp_dim_base, dimension_nonint_exp_neg_base,
    create_save_dir,
    redef_constant, redef_default_constant,
    hashed_none_expr,
    bad_thread_count,
};

std::string unit_to_str(const quantity &q) {
    std::string res = phys::units::to_eng_unit(q);
    if (res.empty()) {
        return "";
    }
    return " " + res;
}


struct dimreal_t {
    mpreal value;
    quantity unit;

    dimreal_t operator+(const dimreal_t &other) const {
        if (!unit.same_dimension(other.unit)) {
            ERR_EXIT(err_t::dimension_add, "attempted to add with different dimension: %s + %s", to_str().c_str(), other.to_str().c_str())
        }
        return dimreal_t{value + other.value, unit};
    }

    dimreal_t operator-(const dimreal_t &other) const {
        if (!unit.same_dimension(other.unit)) {
            ERR_EXIT(err_t::dimension_add, "attempted to subtract with different dimension: %s - %s", to_str().c_str(), other.to_str().c_str())
        }
        return dimreal_t{value - other.value, unit};
    }

    dimreal_t operator-() const {
        return dimreal_t{-value, unit};
    }

    dimreal_t operator*(const dimreal_t &other) const {
        return dimreal_t{value * other.value, unit * other.unit};
    }

    dimreal_t operator/(const dimreal_t &other) const {
        return dimreal_t{value / other.value, unit / other.unit};
    }

    bool operator==(const dimreal_t &other) const {
        return value == other.value && unit == other.unit;
    }

    dimreal_t pow(const dimreal_t &other) const {
        if (other.unit.dimension() != phys::units::dimensionless_d) {
            ERR_EXIT(err_t::dimension_dim_exp, "attempted to exponentiate with non-dimensionless exponent: %s ^ %s", to_str().c_str(), other.to_str().c_str())
        }
        if (!mpfr::isint(other.value) && !unit.same_dimension(quantity())) {
            ERR_EXIT(err_t::dimension_nonint_exp_dim_base, "attempted to exponentiate with non-integer exponent and non-dimensionless base: %s ^ %s", to_str().c_str(), other.to_str().c_str())
        }
        if (unit.dimension() == phys::units::dimensionless_d) {
            if (value < 0 && !mpfr::isint(other.value)) {
                ERR_EXIT(err_t::dimension_nonint_exp_neg_base, "attempted to exponentiate with a non-integer exponent and a negative base: %s ^ %s", to_str().c_str(), other.to_str().c_str())
            }
            return dimreal_t{mpfr::pow(value, other.value), unit}; /* unit is gonna be nothing really here anyways */
        }
        return dimreal_t{mpfr::pow(value, other.value), phys::units::nth_power(unit, static_cast<int>(other.value.toLong()))};
    }

    std::string to_str(int n = digits_prec) const {
        return value.toString(n) + unit_to_str(unit);
    }

    std::string to_str_value(int n = digits_prec) const {
        return value.toString(n);
    }

    static dimreal_t from_str(const std::string &text) {
        auto a = dimreal_t{text, phys::units::to_unit(text, phys::units::dimensionless())};
        std::cout << text << " : " << a.to_str(5) << '\n';
        return a;
    }
};

mpreal cost(const mpreal &a, const mpreal &b) {
    return mpfr::abs(a - b);
}

void split(const std::string &s, const std::string &delim, std::vector<std::string> &outs) {
    std::size_t last = 0, next = 0;
    while ((next = s.find(delim, last)) != std::string::npos) {
        outs.push_back(s.substr(last, next - last));
        last = next + 1;
    }
    if (last != s.size()) { outs.push_back(s.substr(last, s.size())); }
}

std::string get_file_content(const std::string& filename) {
    std::ifstream file(filename);
    std::stringstream ss;
    ss << file.rdbuf();
    std::string text = ss.str();
    file.close();
    return text;
}

struct cnst_t {
    dimreal_t value;
    std::string name;
    bool is_default = false;
};

std::vector<cnst_t> constants;
std::unique_ptr<dimreal_t> target;

enum struct etype_t : std::uint32_t {
    litexpr, cnstexpr,
    addexpr, subexpr,
    mulexpr, divexpr,
    powexpr,
    none
};

struct expr_t { /* NOLINT */
    std::vector<std::shared_ptr<expr_t>> exprs;
    std::vector<std::shared_ptr<expr_t>> parents;
    etype_t type = etype_t::none;
    std::atomic_bool dirty = true;
    dimreal_t cache;

    explicit expr_t() = default;
    explicit expr_t(decltype(exprs) exprs, decltype(parents) parents) : exprs(std::move(exprs)), parents(std::move(parents)) {}
    virtual ~expr_t() = default;

    bool operator==(const expr_t &other) const {
        bool res = true;
        for (std::uint32_t i = 0; i < exprs.size(); i++) {
            res = res && (*exprs[i] == *other.exprs[i]);
        }
        for (std::uint32_t i = 0; i < parents.size(); i++) {
            res = res && (*parents[i] == *other.parents[i]);
        }
        return res;
    }

    std::uint32_t size() {
        std::uint32_t sum = 1;
        for (const std::shared_ptr<expr_t> &expr : exprs) {
            sum += expr->size();
        }
        return sum;
    }

    void signal_dirty() {
        dirty = true;
        for (std::shared_ptr<expr_t> &parent : parents) {
            parent->signal_dirty();
        }
    }

    dimreal_t load() {
        if (dirty) {
            cache = rload();
            dirty = false;
        }
        return cache;
    }
    
    virtual dimreal_t rload() = 0;
    virtual std::string disp() = 0;
};

using sptrexpr_t = std::shared_ptr<expr_t>;

std::vector<std::pair<mpreal, sptrexpr_t>> all;

struct funcexpr_t : virtual expr_t {
    std::string name = "_funcexpr";
    std::function<dimreal_t (const std::vector<sptrexpr_t>&)> f;
    
    explicit funcexpr_t(std::string name, decltype(f) f, std::vector<sptrexpr_t> exprs, std::vector<sptrexpr_t> parents) : expr_t(std::move(exprs), std::move(parents)), name(std::move(name)), f(std::move(f)) {}

    dimreal_t rload() override {
        return f(exprs);
    }

    std::string disp() override {
        std::string arglist;
        for (std::uint32_t i = 0; i < exprs.size(); i++) {
            arglist += exprs[i]->disp() + (i == exprs.size() - 1 ? "" : ", ");
        }
        return name + "(" + arglist + ")";
    }
};

struct unexpr_t : virtual expr_t {
    std::string name = "_unexpr";

    explicit unexpr_t(sptrexpr_t a, std::string name) : name(std::move(name)) {
        exprs.emplace_back(std::move(a));
    }

    
    dimreal_t rload() override {
        return exprs[0]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct litexpr_t : virtual expr_t {
    dimreal_t value;

    explicit litexpr_t(dimreal_t value) : value(std::move(value)) { type = etype_t::litexpr; }

    
    dimreal_t rload() override {
        return value;
    }

    std::string disp() override {
        return value.to_str(digits_prec);
    }
};

struct cnstexpr_t : virtual expr_t {
    std::string name = "_cnstexpr";
    dimreal_t value;

    explicit cnstexpr_t(dimreal_t value, std::string name) : name(std::move(name)), value(std::move(value)) { type = etype_t::cnstexpr; }
    explicit cnstexpr_t(cnst_t constant) : name(std::move(constant.name)), value(constant.value) { type = etype_t::cnstexpr; }

    
    dimreal_t rload() override {
        return value;
    }

    std::string disp() override {
        return name;
    }
};

struct binexpr_t : virtual expr_t {
    std::string name = "_binexpr";

    explicit binexpr_t(sptrexpr_t a, sptrexpr_t b, std::string name) : name(std::move(name)) {
        exprs.emplace_back(std::move(a));
        exprs.emplace_back(std::move(b));
    }
};

struct addexpr_t : public binexpr_t {
    explicit addexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "+") { type = etype_t::addexpr; }

    
    dimreal_t rload() override {
        return exprs[0]->load() + exprs[1]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct subexpr_t : public binexpr_t {
    explicit subexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "-") { type = etype_t::subexpr; }

    
    dimreal_t rload() override {
        return exprs[0]->load() - exprs[1]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct mulexpr_t : public binexpr_t {
    explicit mulexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "*") { type = etype_t::mulexpr; }

    
    dimreal_t rload() override {
        return exprs[0]->load() * exprs[1]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct divexpr_t : public binexpr_t {
    explicit divexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "/") { type = etype_t::divexpr; }

    
    dimreal_t rload() override {
        return exprs[0]->load() / exprs[1]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct powexpr_t : public binexpr_t {
    explicit powexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "^") { type = etype_t::powexpr; }

    
    dimreal_t rload() override {
        return exprs[0]->load().pow(exprs[1]->load());
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

std::uint32_t simplify_passes = 0;

void simplify(sptrexpr_t &a) {
    if (a->size() <= 1) { return; } /* it's just one thing, can't be simplified */

    if (a->type == etype_t::addexpr) {

        /* 0 + expr */
        if (a->exprs[0]->load() == dimreal_t{0}) {
            a = std::move(a->exprs[1]);
            goto redo;
        }

        /* expr + 0 */
        if (a->exprs[1]->load() == dimreal_t{0}) {
            a = std::move(a->exprs[0]);
            goto redo;
        }

    } else if (a->type == etype_t::subexpr) {

        /* expr - 0 */
        if (a->exprs[1]->load() == dimreal_t{0}) {
            a = std::move(a->exprs[0]);
            goto redo;
        }

    } else if (a->type == etype_t::mulexpr) {

        /* 1 * expr */
        if (a->exprs[0]->load() == dimreal_t{1}) {
            a = std::move(a->exprs[1]);
            goto redo;
        }

        /* expr * 1 */
        if (a->exprs[1]->load() == dimreal_t{1}) {
            a = std::move(a->exprs[0]);
            goto redo;
        }

    } else if (a->type == etype_t::divexpr) {

        /* expr / 1 */
        if (a->exprs[1]->load() == dimreal_t{1}) {
            a = std::move(a->exprs[0]);
            goto redo;
        }

        /* 1 / expr */
        if (a->exprs[0]->load() == dimreal_t{1}) {

            /* 1 / (expr / expr) */
            if (a->exprs[1]->type == etype_t::divexpr) { /* invert it */
                sptrexpr_t tmp = std::move(a->exprs[1]->exprs[0]);
                a->exprs[1]->exprs[0] = std::move(a->exprs[1]->exprs[1]);
                a->exprs[1]->exprs[1] = std::move(tmp);
                a = std::move(a->exprs[1]);
                goto redo;
            }

        }

    }

    for (sptrexpr_t &expr : a->exprs) {
        simplify(expr);
    }
    return;

    redo: {
        simplify(a);
    }
}


void recurse(sptrexpr_t, std::uint32_t);

void test_expr(const sptrexpr_t& a, std::uint32_t cursize) {
    const dimreal_t res = a->load();
    const mpreal diff = cost(res.value, target->value);
    if (res.unit.same_dimension(target->unit)) {
        all.emplace_back(diff, a);
    }
    recurse(a, cursize + 1);
}


void recurse(sptrexpr_t b, std::uint32_t cursize = 1) {
    if (cursize > max_expr_size) { return; }
    for (const cnst_t &constant : constants) {
        /* don't technically need to include constant - b or constant / b, 
         * it's covered by the 0 - and 1 / cases in the next recursion
         * however, we're not guaranteed another recursion due to limits on expr size
         * so we do them here anyway
         */
        if (b->load().unit.same_dimension(quantity())) {
            if (constant.value.unit.same_dimension(quantity())) {
                if (constant.value.value > 0 || (constant.value.value < 0 && mpfr::isint(b->load().value))) {
                    test_expr(std::make_shared<powexpr_t>(std::make_shared<cnstexpr_t>(constant), b), cursize);
                }
                if (b->load().value > 0 || (b->load().value < 0 && mpfr::isint(constant.value.value))) {
                    test_expr(std::make_shared<powexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
                }
                if (mpfr::isint(constant.value.value)) {
                    test_expr(std::make_shared<powexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
                }
            }
            if (mpfr::isint(b->load().value)) {
                test_expr(std::make_shared<powexpr_t>(std::make_shared<cnstexpr_t>(constant), b), cursize);
            }
        }
        test_expr(std::make_shared<mulexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
        if (constant.value.value != 0) {
            test_expr(std::make_shared<divexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
        }
        if (b->load().value != 0) {
            test_expr(std::make_shared<divexpr_t>(std::make_shared<cnstexpr_t>(constant), b), cursize);
        }
        if (constant.value.unit.same_dimension(b->load().unit)) {
            test_expr(std::make_shared<addexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
            test_expr(std::make_shared<subexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
            test_expr(std::make_shared<subexpr_t>(std::make_shared<cnstexpr_t>(constant), b), cursize);
        }
    }

    for (mpreal i = 2; i <= max_int_constants; i++) {
        test_expr(std::make_shared<mulexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i, target->unit / b->load().unit})), cursize);
        sptrexpr_t bottom = std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit / target->unit});
        if (bottom->load().value != 0) {
            test_expr(std::make_shared<divexpr_t>(b, bottom), cursize);
        }
        if (b->load().unit.same_dimension(quantity())) {
            if (mpfr::isint(b->load().value)) {
                test_expr(std::make_shared<powexpr_t>(std::make_shared<litexpr_t>(dimreal_t{i}), b), cursize);
            }
            test_expr(std::make_shared<powexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i})), cursize);
        }
    }

    for (mpreal i = 1; i <= max_int_constants; i++) {
        if (b->load().value != 0) {
            test_expr(std::make_shared<divexpr_t>(std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit * target->unit}), b), cursize);
        }
        test_expr(std::make_shared<addexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit})), cursize);
        test_expr(std::make_shared<subexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit})), cursize);
        test_expr(std::make_shared<subexpr_t>(std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit}), b), cursize);
    }
    test_expr(std::make_shared<subexpr_t>(std::make_shared<litexpr_t>(dimreal_t{0, b->load().unit}), b), cursize);
}


int main(int argc, char **argv) {
    std::int32_t thread_count = 1;

    if (argc == 2) {
        if (!std::strcmp(argv[1], "-v") || !std::strcmp(argv[1], "--version")) {
            std::cout << "exactonator version " VERSION ", " YEAR " by .stole.\n";
            return 0;
        }

        if (!std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) {
            std::cout << "usage: " << argv[0] << R"( [flags]

    -j <count> : runs <count> threads
    -v, --version : displays texproj's version
    -h, --help : displays this help
)";
            return 0;
        }
    } else if (argc > 2) {
        for (std::uint32_t i = 1; i < argc - 1; i++) {
            if (argv[i][0] != '-') {
                std::cerr << "error: unexpected option \"" << argv[i] << "\"\n";
                return 4;
            }
            if (!std::strcmp(argv[i], "-j")) {
                thread_count = std::strtol(argv[++i], nullptr, 0);
                if (thread_count == 0) {
                    ERR_EXIT(err_t::bad_thread_count, "bad thread count, must be an integer > 0")
                }
            }
        }
    }

    if (!std::filesystem::exists(SAVE_AST_DIR)) {
        if (!std::filesystem::create_directory(SAVE_AST_DIR)) {
            ERR_EXIT(err_t::create_save_dir, "AST save directory \"%s\" does not exist, failed to create", SAVE_AST_DIR)
        }
    }

    std::cout << "digits: ";
    std::string digits_prec_str;
    std::getline(std::cin, digits_prec_str);
    digits_prec = std::stoi(digits_prec_str);
    mpreal::set_default_prec(static_cast<mpfr_prec_t>(std::pow(10, std::log2(digits_prec + 1))));

    std::cout << "target: ";
    std::string target_str;
    std::getline(std::cin, target_str);
    target = std::make_unique<dimreal_t>(dimreal_t::from_str(target_str));

    std::cout << "max expr size: ";
    std::string max_expr_size_str;
    std::getline(std::cin, max_expr_size_str);
    max_expr_size = std::stoi(max_expr_size_str);

    std::cout << "integer constants up to: ";
    std::string max_int_constants_str;
    std::getline(std::cin, max_int_constants_str);
    max_int_constants = std::stoi(max_int_constants_str);

    const std::vector<cnst_t> default_constants = {
        cnst_t{
            dimreal_t{mpfr::const_pi()},
            "pi", true
        },
        cnst_t{
            dimreal_t{const_e_str},
            "e", true
        },
        cnst_t{
            dimreal_t{mpfr::const_euler()},
            "euler", true
        },
        cnst_t{
            dimreal_t{mpfr::const_log2()},
            "ln2", true
        },
        cnst_t{
            dimreal_t{mpfr::const_catalan()},
            "catalan", true
        },
        cnst_t{
            dimreal_t{("1" + mpfr::sqrt("5")) / "2"},
            "phi", true
        },
        cnst_t{
            dimreal_t{"0.0072973525693"},
            "fine-structure", true
        }
    };

    std::ifstream constants_file(CONSTANTS_FILENAME);
    std::uint32_t ocount = 0;
    for (std::string line; std::getline(constants_file, line);) {
        if (line.empty()) { continue; }
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        ocount++;
        std::vector<std::string> opts;
        opts.reserve(2);
        split(line, "=", opts);

        if (opts.size() == 1) {
            auto pos = std::find_if(default_constants.begin(), default_constants.end(), [&](const cnst_t &constant) -> bool { return constant.name == opts[0]; });
            if (pos != default_constants.end()) {
                constants.push_back(*pos);
                continue;
            }
            std::cout << "warning: " << CONSTANTS_FILENAME << " file #" << ocount << " config option had 1 token; expecting default constant name, but \"" << opts[0] << "\" is not of {";
            for (std::uint32_t i = 0; i < default_constants.size() - 1; i++) {
                std::cout << "\"" << default_constants[i].name << "\", ";
            }
            std::cout << "\"" << (default_constants.end() - 1)->name << "\"}. specifying a value might look like \"" << opts[0] << " = " << "1.0 s\" ... skipping\n";
            continue;
        } else if (opts.size() < 2) {
            std::cout << "warning: " << CONSTANTS_FILENAME << " file #" << ocount << " config option had " << opts.size() << " token" << (opts.size() == 1 ? "" : "s") << " ... skipping\n";
            continue;
        }

        if (opts.size() > 2) {
            std::cout << "warning: " << CONSTANTS_FILENAME << " file #" << ocount << " config option had " << opts.size() << " tokens ... using first two\n";
        }

        std::string name = opts[0], value = opts[1];
        if (name.empty() || value.empty()) {
            std::cout << "warning: " << CONSTANTS_FILENAME << " file #" << ocount << " config option had empty name or value ... skipping\n";
            continue;
        }
        for (std::uint32_t i = 0; i < constants.size(); i++) {
            if (constants[i].name == name) {
                ERR_EXIT(err_t::redef_constant, "\"%s\" file #%i config option redefined: %s = %s over previous #%i config option %s = %s\n", CONSTANTS_FILENAME, ocount, name.c_str(), dimreal_t::from_str(value).to_str().c_str(), i + 1, constants[i].name.c_str(), constants[i].value.to_str().c_str())
            }
        }
        for (std::uint32_t i = 0; i < default_constants.size(); i++) {
            if (default_constants[i].name == name) {
                ERR_EXIT(err_t::redef_default_constant, "\"%s\" file #%i config option redefined default constant: %s = %s over %s = %s\n", CONSTANTS_FILENAME, ocount, name.c_str(), dimreal_t::from_str(value).to_str().c_str(), default_constants[i].name.c_str(), constants[i].value.to_str().c_str())
            }
        }
        /* std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::tolower(c); }); */

        cnst_t tcnst = cnst_t{.value = dimreal_t::from_str(value), .name = name};

        constants.push_back(tcnst);
    }


    /* this section shouldn't be changed */
    /* ---- */
    std::stringstream seed_str_stream;
    seed_str_stream << "max_expr=" << max_expr_size << ",max_int=" << max_int_constants << ';';
    if (!constants.empty()) {
        std::uint32_t default_constant_count = 0;
        for (std::uint32_t i = 0; i < constants.size() - 1; i++) {
            const cnst_t &constant = constants[i];
            if (constant.is_default) {
                seed_str_stream << constant.name;
            } else {
                seed_str_stream << "%" << default_constant_count++ << "=" << constant.value.to_str();
            }
            seed_str_stream << ",";
        }
        const cnst_t &constant = *(constants.end() - 1);
        if (constant.is_default) {
            seed_str_stream << constant.name;
        } else {
            seed_str_stream << "%" << constants.size() - 1 << "=" << constant.value.to_str();
        }
    }

    std::string seed_str = seed_str_stream.str();

    std::stringstream savefilename;
    savefilename << std::hex << std::hash<std::string>{}(seed_str);

    std::filesystem::current_path(SAVE_AST_DIR);
    std::ofstream savefile(savefilename.str());
    savefile << seed_str << '\n';
    savefile.close();
    /* ---- */

    for (const cnst_t &constant : constants) {
        test_expr(std::make_shared<cnstexpr_t>(constant), 1);
    }

    for (mpreal i = 1; i <= max_int_constants; i++) {
        test_expr(std::make_shared<litexpr_t>(dimreal_t{i, target->unit}), 1);
    }

    decltype(all) selected;

    /* remove expressions with the same error, picking the one with the smallest expression length */
    for (const std::pair<mpreal, sptrexpr_t> &a : all) {
        auto pos = std::find_if(selected.begin(), selected.end(), [&](const std::pair<mpreal, sptrexpr_t> &b) { return a.first == b.first; });
        if (pos != selected.end()) {
            if (pos->second->size() > a.second->size()) {
                *pos = a;
            }
        } else {
            selected.emplace_back(a.first, a.second);
        }
    }

    std::sort(selected.begin(), selected.end(), [](const std::pair<mpreal, sptrexpr_t> &a, const std::pair<mpreal, sptrexpr_t> &b) {
        return a.first < b.first;
    });

    for (std::uint32_t i = 0; i < 30 && i < selected.size(); i++) {
        std::cout << selected[i].second->disp() << " | err: " << selected[i].first.toString(digits_prec) << '\n';
    }

    return 0;
}
