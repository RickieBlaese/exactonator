#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <unordered_map>
#include <utility>
#include <optional>
#include <vector>
#include <stdexcept>

#include <cinttypes>

#include "mpreal/mpreal.h"

#include "const_e.h"


#include "phys/units/io.hpp"
#include "phys/units/quantity.hpp"
#include "phys/units/other_units.hpp"

using mpfr::mpreal;
using phys::units::quantity;

std::int32_t digits_prec = 0;
std::int32_t max_int_constants = 0;
std::int32_t max_expr_size = 1;

#define ERR_EXIT(...) { \
    std::fprintf(stderr, "error: file " __FILE__ ":%i in %s(): ", __LINE__, __func__); \
    std::fprintf(stderr, __VA_ARGS__); \
    std::fputc('\n', stderr); \
    throw std::runtime_error(""); \
}

static constexpr const char CONSTANTS_FILENAME[] = "constants.conf";

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
            ERR_EXIT("attempted to add with different dimension: %s + %s", to_str().c_str(), other.to_str().c_str());
        }
        return dimreal_t{value + other.value, unit};
    }

    dimreal_t operator-(const dimreal_t &other) const {
        if (!unit.same_dimension(other.unit)) {
            ERR_EXIT("attempted to subtract with different dimension: %s - %s", to_str().c_str(), other.to_str().c_str());
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
            ERR_EXIT("attempted to exponentiate with non-dimensionless exponent: %s ^ %s", to_str().c_str(), other.to_str().c_str());
        }
        if (!mpfr::isint(other.value) && !unit.same_dimension(quantity())) {
            ERR_EXIT("attempted to exponentiate with non-integer exponent and non-dimensionless base: %s ^ %s", to_str().c_str(), other.to_str().c_str());
        }
        if (unit.dimension() == phys::units::dimensionless_d) {
            if (value < 0 && !mpfr::isint(other.value)) {
                ERR_EXIT("attempted to exponentiate with a non-integer exponent and a negative base: %s ^ %s", to_str().c_str(), other.to_str().c_str());
            }
            return dimreal_t{mpfr::pow(value, other.value), unit}; /* unit is gonna be nothin rly here anyways */
        }
        return dimreal_t{mpfr::pow(value, other.value), phys::units::nth_power(unit, static_cast<int>(other.value.toLong()))};
    }

    std::string to_str(int n = digits_prec) const {
        return value.toString(n) + unit_to_str(unit);
    }

    static dimreal_t from_str(const std::string &text) {
        return dimreal_t{text, phys::units::to_unit(text, phys::units::dimensionless())};
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
};

std::vector<cnst_t> constants;
std::optional<dimreal_t> target;


enum struct etype_t : std::uint32_t {
    litexpr, cnstexpr,
    addexpr, subexpr,
    mulexpr, divexpr,
    powexpr,
    none
};

struct expr_t {
    std::vector<std::shared_ptr<expr_t>> exprs;
    etype_t type = etype_t::none;

    explicit expr_t() = default;
    explicit expr_t(decltype(exprs) exprs) : exprs(std::move(exprs)) {}
    virtual ~expr_t() = default;

    bool operator==(const expr_t &other) const {
        bool res = true;
        for (std::uint32_t i = 0; i < exprs.size(); i++) {
            res = res && (*exprs[i] == *other.exprs[i]);
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
    
    virtual dimreal_t load() = 0;
    virtual quantity get_unit() = 0;
    virtual std::string disp() = 0;
};

using sptrexpr_t = std::shared_ptr<expr_t>;

std::vector<std::pair<mpreal, sptrexpr_t>> all;

struct unexpr_t : virtual expr_t {
    std::string name = "_unexpr";

    explicit unexpr_t(sptrexpr_t a, std::string name) : name(std::move(name)) {
        exprs.emplace_back(std::move(a));
    }
};

struct litexpr_t : virtual expr_t {
    dimreal_t value;

    explicit litexpr_t(dimreal_t value) : value(std::move(value)) { type = etype_t::litexpr; }

    dimreal_t load() override {
        return value;
    }

    quantity get_unit() override {
        return value.unit;
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

    dimreal_t load() override {
        return value;
    }

    quantity get_unit() override {
        return value.unit;
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

    quantity get_unit() override {
        return exprs[0]->get_unit(); /* if they're different, then it is malformed anyways, and we'll hit an error later */
    }

};

struct addexpr_t : public binexpr_t {
    explicit addexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "+") { type = etype_t::addexpr; }

    dimreal_t load() override {
        return exprs[0]->load() + exprs[1]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct subexpr_t : public binexpr_t {
    explicit subexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "-") { type = etype_t::subexpr; }

    dimreal_t load() override {
        return exprs[0]->load() - exprs[1]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct mulexpr_t : public binexpr_t {
    explicit mulexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "*") { type = etype_t::mulexpr; }

    dimreal_t load() override {
        return exprs[0]->load() * exprs[1]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct divexpr_t : public binexpr_t {
    explicit divexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "/") { type = etype_t::divexpr; }

    dimreal_t load() override {
        return exprs[0]->load() / exprs[1]->load();
    }

    std::string disp() override {
        return "(" + exprs[0]->disp() + " " + name + " " + exprs[1]->disp() + ")";
    }
};

struct powexpr_t : public binexpr_t {
    explicit powexpr_t(sptrexpr_t a, sptrexpr_t b) : binexpr_t(std::move(a), std::move(b), "^") { type = etype_t::powexpr; }

    dimreal_t load() override {
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
    if (a->load().unit.same_dimension(target.value().unit)) {
        all.emplace_back(cost(a->load().value, target.value().value), a);
    }
    recurse(a, cursize + 1);
}


void recurse(sptrexpr_t b, std::uint32_t cursize = 1) {
    if (cursize > max_expr_size) { return; }
    for (const cnst_t &constant : constants) {
        /* don't technically need to include constant - b or constant / b, 
         * it's covered by the 0 - and 1 / cases in the next recursion.
         * however, we're not guaranteed another recursion due to limits on expr size.
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
            }
            if (mpfr::isint(b->load().value)) {
                test_expr(std::make_shared<powexpr_t>(std::make_shared<cnstexpr_t>(constant), b), cursize);
            }
            if (mpfr::isint(constant.value.value)) {
                test_expr(std::make_shared<powexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
            }
        }
        if (constant.value.unit.same_dimension(b->load().unit)) {
            test_expr(std::make_shared<addexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
            test_expr(std::make_shared<subexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
            test_expr(std::make_shared<subexpr_t>(std::make_shared<cnstexpr_t>(constant), b), cursize);
        }
        test_expr(std::make_shared<mulexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
        test_expr(std::make_shared<divexpr_t>(b, std::make_shared<cnstexpr_t>(constant)), cursize);
        test_expr(std::make_shared<divexpr_t>(std::make_shared<cnstexpr_t>(constant), b), cursize);
    }

    for (mpreal i = 2; i < max_int_constants; i++) {
        test_expr(std::make_shared<mulexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i, target.value().unit / b->load().unit})), cursize);
        test_expr(std::make_shared<divexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit / target.value().unit})), cursize);
        if (b->load().unit.same_dimension(quantity())) {
            if (mpfr::isint(b->load().value)) {
                test_expr(std::make_shared<powexpr_t>(std::make_shared<litexpr_t>(dimreal_t{i}), b), cursize);
            }
            test_expr(std::make_shared<powexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i})), cursize);
        }
    }

    for (mpreal i = 1; i < max_int_constants; i++) {
        test_expr(std::make_shared<addexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit})), cursize);
        test_expr(std::make_shared<subexpr_t>(b, std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit})), cursize);
        test_expr(std::make_shared<subexpr_t>(std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit}), b), cursize);
        test_expr(std::make_shared<divexpr_t>(std::make_shared<litexpr_t>(dimreal_t{i, b->load().unit * target.value().unit}), b), cursize);
    }
}


int main() {
    std::cout << "digits: ";
    std::cin >> digits_prec;
    mpreal::set_default_prec(static_cast<mpfr_prec_t>(std::pow(10, std::log2(digits_prec + 1))));

    std::cout << "target: ";
    std::string target_str;
    std::cin >> target_str;
    target = dimreal_t::from_str(target_str);

    std::cout << "max expr size: ";
    std::cin >> max_expr_size;

    std::cout << "integer constants up to: ";
    std::cin >> max_int_constants;

    std::ifstream constants_file(CONSTANTS_FILENAME);
    std::uint32_t ocount = 0;
    for (std::string line; std::getline(constants_file, line);) {
        if (line.empty()) { continue; }
        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
        ocount++;
        std::vector<std::string> opts;
        opts.reserve(2);
        split(line, "=", opts);

        if (opts.size() < 2) {
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
                std::cerr << "fatal: " << CONSTANTS_FILENAME << " file #" << ocount << " config option redefined: " << name << " = " << value << " over previous #" << i + 1 << " config option " << constants[i].name << " = " << constants[i].value.to_str() << '\n';
                return 1;
            }
        }
        /* std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::tolower(c); }); */

        cnst_t tcnst = cnst_t{.value = dimreal_t::from_str(value), .name = name};

        constants.push_back(tcnst);
    }

    for (const cnst_t &constant : constants) {
        test_expr(std::make_shared<cnstexpr_t>(constant), 1);
    }

    for (mpreal i = 1; i < max_int_constants; i++) {
        test_expr(std::make_shared<litexpr_t>(dimreal_t{i, target.value().unit}), 1);
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
        std::cout << selected[i].second->disp() << " | err: " << selected[i].first << '\n';
    }


    return 0;
}
