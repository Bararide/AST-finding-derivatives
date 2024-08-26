#include <functional>
#include <complex>
#include <tuple>
#include <string>
#include <type_traits>
#include <memory>
#include <stack>
#include <vector>
#include <map>
#include <typeinfo>
#include <iostream>

using value_t = std::complex<double>;
using func_t = std::function<value_t(value_t)>;

namespace AST {
    class ASTNode {
    public:
        virtual ~ASTNode() = default;
        virtual std::string toString() const = 0;
        
        virtual long double evaluate(const std::map<std::string, long double>& variables) const = 0;
        virtual value_t evaluateComplex(const std::map<std::string, value_t>& variables) const = 0;

        virtual std::shared_ptr<ASTNode> derivative(const std::string& var) const = 0;
        virtual std::shared_ptr<ASTNode> clone() const = 0;
    };

    class NumberNode : public ASTNode {
        long double value;
    public:
        NumberNode(long double val) : value(val) {}
        std::string toString() const override {
            return std::to_string(value);
        }
        long double evaluate(const std::map<std::string, long double>& variables) const override {
            return value;
        }
        value_t evaluateComplex(const std::map<std::string, value_t>& variables) const override {
            return value_t(value, 0.L);
        }

        std::shared_ptr<ASTNode> derivative(const std::string& var) const override {
            return std::make_shared<NumberNode>(0.L);
        }

        std::shared_ptr<ASTNode> clone() const override {
            return std::make_shared<NumberNode>(value);
        }
    };

    class VariableNode : public ASTNode {
        std::string name;
    public:
        VariableNode(const std::string& n) : name(n) {}
        std::string toString() const override {
            return name;
        }

        long double evaluate(const std::map<std::string, long double>& variables) const override {
            auto it = variables.find(name);
            if (it == variables.end()) {
                throw std::runtime_error("Undefined variable: " + name);
            }
            return it->second;
        }

        value_t evaluateComplex(const std::map<std::string, value_t>& variables) const override {
            auto it = variables.find(name);
            if (it == variables.end()) {
                throw std::runtime_error("Undefined variable: " + name);
            }
            return it->second;
        }

        std::shared_ptr<ASTNode> derivative(const std::string& var) const override {
            if (name == var) {
                return std::make_shared<NumberNode>(1.L);
            } else {
                return std::make_shared<NumberNode>(0.L);
            }
        }

        std::shared_ptr<ASTNode> clone() const override {
            return std::make_shared<VariableNode>(name);
        }
    };

    class OpNode : public ASTNode {
        std::string op;
        std::shared_ptr<ASTNode> operand1;
        std::shared_ptr<ASTNode> operand2;

    public:
        OpNode(const std::string& o, std::shared_ptr<ASTNode> node1)
            : op(o), operand1(std::move(node1)), operand2(nullptr) {}

        OpNode(const std::string& o, std::shared_ptr<ASTNode> node1, std::shared_ptr<ASTNode> node2)
            : op(o), operand1(std::move(node1)), operand2(std::move(node2)) {}

        std::string toString() const override {
            if (operand2) {
                return "(" + operand1->toString() + " " + op + " " + operand2->toString() + ")";
            } else {
                return op + "(" + operand1->toString() + ")";
            }
        }

        long double evaluate(const std::map<std::string, long double>& variables) const override {
            long double val1 = operand1->evaluate(variables);
            if (operand2) {
                long double val2 = operand2->evaluate(variables);
                if (op == "+") return val1 + val2;
                if (op == "-") return val1 - val2;
                if (op == "*") return val1 * val2;
                if (op == "/") return val1 / val2;
                if (op == "^") return std::pow(val1, val2);
                throw std::runtime_error("Unknown binary operator: " + op);
            } else {
                if (op == "sin")  return std::sin(val1);
                if (op == "cos")  return std::cos(val1);
                if (op == "tan")  return std::tan(val1);
                if (op == "sqrt") return std::sqrt(val1);
                if (op == "log")  return std::log(val1);
                if (op == "cot")  return 1.0L / std::tan(val1);
                if (op == "ln")   return std::log1pl(val1);
                throw std::runtime_error("Unknown unary operator: " + op);
            }
        }

        value_t evaluateComplex(const std::map<std::string, value_t>& variables) const override {
            value_t val1 = operand1->evaluateComplex(variables);
            if (operand2) {
                value_t val2 = operand2->evaluateComplex(variables);
                if (op == "+") return val1 + val2;
                if (op == "-") return val1 - val2;
                if (op == "*") return val1 * val2;
                if (op == "/") return val1 / val2;
                if (op == "^") return std::pow(val1, val2);
                throw std::runtime_error("Unknown binary operator: " + op);
            } else {
                if (op == "sin")  return std::sin(val1);
                if (op == "cos")  return std::cos(val1);
                if (op == "tan")  return std::tan(val1);
                if (op == "sqrt") return std::sqrt(val1);
                if (op == "log")  return std::log(val1);
                if (op == "cot")  return 1.0 / std::tan(val1);
                throw std::runtime_error("Unknown unary operator: " + op);
            }
        }

        std::shared_ptr<ASTNode> clone() const override {
            if (operand2) {
                return std::make_shared<OpNode>(op, operand1->clone(), operand2->clone());
            }
            else {
                return std::make_shared<OpNode>(op, operand1->clone());
            }
        }

        std::shared_ptr<ASTNode> derivative(const std::string& var) const override {
            if (operand2) {
                // Для бинарных операций
                if (op == "+") {
                    return std::make_shared<OpNode>("+", operand1->derivative(var), operand2->derivative(var));
                }
                else if (op == "-") {
                    return std::make_shared<OpNode>("-", operand1->derivative(var), operand2->derivative(var));
                }
                else if (op == "*") {
                    return std::make_shared<OpNode>("+",
                        std::make_shared<OpNode>("*", operand1->clone(), operand2->derivative(var)),
                        std::make_shared<OpNode>("*", operand1->derivative(var), operand2->clone())
                    );
                }
                else if (op == "/") {
                    return std::make_shared<OpNode>("/",
                        std::make_shared<OpNode>("-", 
                            std::make_shared<OpNode>("*", operand2->clone(), operand1->derivative(var)),
                            std::make_shared<OpNode>("*", operand1->clone(), operand2->derivative(var))),
                        std::make_shared<OpNode>("^", operand2->clone(), std::make_shared<NumberNode>(2))
                    );
                }
                else if (op == "^") {
                    // d/dx(u^v) = u^v * (v * du/dx / u + ln(u) * dv/dx)
                    return std::make_shared<OpNode>("*",
                        std::make_shared<OpNode>("^", operand1->clone(), operand2->clone()),
                        std::make_shared<OpNode>("+",
                            std::make_shared<OpNode>("*",
                                std::make_shared<OpNode>("/",
                                    std::make_shared<OpNode>("*", operand2->clone(), operand1->derivative(var)),
                                    operand1->clone()
                                ),
                                std::make_shared<NumberNode>(1.0)
                            ),
                            std::make_shared<OpNode>("*",
                                std::make_shared<OpNode>("log", operand1->clone()),
                                operand2->derivative(var)
                            )
                        )
                    );
                }
            }
            else {
                if (op == "sin") {
                    return std::make_shared<OpNode>("*", std::make_shared<OpNode>("cos", operand1->clone()), operand1->derivative(var));
                }
                else if (op == "cos") {
                    return std::make_shared<OpNode>("*", std::make_shared<OpNode>("*",
                        std::make_shared<NumberNode>(-1),
                        std::make_shared<OpNode>("sin", operand1->clone())),
                        operand1->derivative(var));
                }
                else if (op == "sqrt") {
                    return std::make_shared<OpNode>("*",
                        std::make_shared<OpNode>("/", std::make_shared<NumberNode>(1),
                        std::make_shared<OpNode>("*", std::make_shared<NumberNode>(2),
                        std::make_shared<OpNode>("sqrt", operand1->clone()))),
                        operand1->derivative(var));
                }
                else if (op == "tan") {
                    // d/dx(tan(u)) = (1 / cos^2(u)) * du/dx
                    return std::make_shared<OpNode>("*",
                        std::make_shared<OpNode>("/",
                            std::make_shared<NumberNode>(1),
                            std::make_shared<OpNode>("^",
                                std::make_shared<OpNode>("cos", operand1->clone()),
                                std::make_shared<NumberNode>(2))),
                        operand1->derivative(var)); 
                }
                else if (op == "cot") {
                    // d/dx(cot(u)) = -(1 / sin^2(u)) * du/dx
                    return std::make_shared<OpNode>("*",
                        std::make_shared<OpNode>("*",
                            std::make_shared<NumberNode>(-1),
                            std::make_shared<OpNode>("/",
                                std::make_shared<NumberNode>(1),
                                std::make_shared<OpNode>("^",
                                    std::make_shared<OpNode>("sin", operand1->clone()),
                                    std::make_shared<NumberNode>(2)))),
                        operand1->derivative(var)); 
                }
                else if (op == "log") {
                    // d/dx(log(u)) = (1/u) * du/dx
                    return std::make_shared<OpNode>("*",
                        std::make_shared<OpNode>("/",
                            std::make_shared<NumberNode>(1),
                            operand1->clone()), // u
                        operand1->derivative(var)); // du/dx
                }
            }

            throw std::runtime_error("Unknown operator for derivation: " + op);
        }
    };

    class ExpressionParser {
    private:
        std::vector<std::string> expr;
        std::string digit = "";

        std::map<std::string, void (*)(std::vector<std::string>&, std::string&)> commands = {
            {"sin", [](std::vector<std::string>& polish, std::string& command) {
                polish.push_back(std::move(command));
            }},
            {"cos", [](std::vector<std::string>& polish, std::string& command) {
                polish.push_back(std::move(command));
            }},
            {"tan", [](std::vector<std::string>& polish, std::string& command) {
                polish.push_back(std::move(command));
            }},
            {"cot", [](std::vector<std::string>& polish, std::string& command) {
                polish.push_back(std::move(command));
            }},
            {"sqrt", [](std::vector<std::string>& polish, std::string& command) {
                polish.push_back(std::move(command));
            }},
            {"log", [](std::vector<std::string>& polish, std::string& command) {
                polish.push_back(std::move(command));
            }}
        };

    public:
        ExpressionParser() {}

        std::shared_ptr<ASTNode> buildAST() {
            std::vector<std::shared_ptr<ASTNode>> stack;
            std::string invalidToken;

            for (auto v : expr) {
                invalidToken += v;
                invalidToken += " ";
            }

            for (const auto& token : expr) {
                try {
                    if (token == "+" || token == "-" || token == "*" || token == "/" || token == "^") {
                        if (stack.size() < 2) {
                            throw std::runtime_error("Invalid expression: " + token);
                        }
                        auto right = std::move(stack.back());
                        stack.pop_back();
                        auto left = std::move(stack.back());
                        stack.pop_back();
                        stack.push_back(std::make_shared<OpNode>(token, std::move(left), std::move(right)));
                    } 
                    else if (token == "sin" || token == "cos" || token == "tan" || token == "cot" || token == "sqrt" || token == "log") {
                        if (stack.empty()) {
                            throw std::runtime_error("Invalid expression: " + token);
                        }
                        auto operand = std::move(stack.back());
                        stack.pop_back();
                        stack.push_back(std::make_shared<OpNode>(token, std::move(operand)));
                    } 
                    else {
                        try {
                            long double value = std::stod(token);
                            stack.push_back(std::make_shared<NumberNode>(value));
                        } 
                        catch (const std::invalid_argument&) {
                            stack.push_back(std::make_shared<VariableNode>(token));
                        }
                    }
                } 
                catch (const std::runtime_error& e) {
                    std::cout << invalidToken << std::endl;
                    throw;
                }
            }

            if (stack.size() != 1) {
                throw std::runtime_error("Invalid expression: " + invalidToken);
            }

            return std::move(stack.back());
        }

        long double evaluateExpression(const std::shared_ptr<ASTNode>& ast, const std::map<std::string, long double>& variables) {
            return ast->evaluate(variables);
        }

        value_t evaluateComplexExpression(const std::shared_ptr<ASTNode>& ast, const std::map<std::string, value_t>& variables) {
            return ast->evaluateComplex(variables);
        }

        std::shared_ptr<ASTNode> computeDerivative(const std::shared_ptr<ASTNode>& ast, const std::string& var) {
            return ast->derivative(var);
        }

        void to_polish_notation(const std::string& expression) {
            std::map<char, std::pair<int, bool>> precedence {
                {'+', {1, false}}, {'-', {1, false}},
                {'*', {2, false}}, {'/', {2, false}},
                {'^', {4, true}}
            };
            std::stack<std::string> oper_stack;
            std::string token = "";

            auto push_token = [&]() {
                if (!token.empty()) {
                    expr.push_back(token);
                    token = "";
                }
            };

            for (size_t i = 0; i < expression.size(); ++i) {
                char c = expression[i];

                if (std::isdigit(c) || c == '.' || (c == '-' && (i == 0 || expression[i-1] == '('))) {
                    token += c;
                    if (i + 1 < expression.size()) {
                        char next = expression[i + 1];
                        if (!std::isdigit(next) && next != '.' && !(c == '-' && std::isdigit(next))) {
                            push_token();
                        }
                    } else {
                        push_token();
                    }
                }
                else if (std::isalpha(c)) {
                    token += c;
                    if (i + 1 == expression.size() || !std::isalpha(expression[i + 1])) {
                        if (commands.find(token) != commands.end()) {
                            oper_stack.push(token);
                        } else {
                            expr.push_back(token);
                        }
                        token = "";
                    }
                }
                else if (c == '(') {
                    oper_stack.push(std::string(1, c));
                }
                else if (c == ')') {
                    push_token();
                    while (!oper_stack.empty() && oper_stack.top() != "(") {
                        expr.push_back(oper_stack.top());
                        oper_stack.pop();
                    }
                    if (!oper_stack.empty() && oper_stack.top() == "(") {
                        oper_stack.pop();
                    }
                    if (!oper_stack.empty() && commands.find(oper_stack.top()) != commands.end()) {
                        expr.push_back(oper_stack.top());
                        oper_stack.pop();
                    }
                }
                else if (precedence.find(c) != precedence.end()) {
                    push_token();
                    while (!oper_stack.empty() && oper_stack.top().size() == 1 && 
                            (precedence[oper_stack.top()[0]].first > precedence[c].first ||
                            (precedence[oper_stack.top()[0]].first == precedence[c].first && 
                            !precedence[c].second))) {
                        expr.push_back(oper_stack.top());
                        oper_stack.pop();
                    }
                    oper_stack.push(std::string(1, c));
                }
            }

            push_token();

            while (!oper_stack.empty()) {
                expr.push_back(oper_stack.top());
                oper_stack.pop();
            }
        }
    };
}

std::tuple<func_t, func_t, func_t> differentiate(const std::string& eq) {
    using namespace AST;

    auto ex = std::make_shared<ExpressionParser>();

    ex->to_polish_notation(eq);
    auto ast = ex->buildAST();

    auto derivative = ex->computeDerivative(ast, "x");
    auto derivative_2 = ex->computeDerivative(derivative, "x");

    return {
        [ex, ast = std::move(ast)](value_t x) {
            std::map<std::string, value_t> variables{{"x", x}};
            return ex->evaluateComplexExpression(ast, variables);
        },
        [ex, derivative = std::move(derivative)](value_t x) {
            std::map<std::string, value_t> variables{{"x", x}};
            return ex->evaluateComplexExpression(derivative, variables);
        },
        [ex, derivative_2 = std::move(derivative_2)](value_t x) {
            std::map<std::string, value_t> variables{{"x", x}};
            return ex->evaluateComplexExpression(derivative_2, variables);
        }
    };
}

int main() {
    using namespace AST;

    ExpressionParser* ex = new ExpressionParser();

    ex->to_polish_notation("4 * log(x) + x^2 / 2^x");

    auto ast = ex->buildAST();

    auto derivative = ex->computeDerivative(ast, "x");

    auto derivative_2 = ex->computeDerivative(derivative, "x");

    const auto [f_x, df_x, df2_x] = differentiate("4 * log(x) + x^2 / 2^x");

    std::cout << f_x({1, 1}) << std::endl;

    return 0;
}