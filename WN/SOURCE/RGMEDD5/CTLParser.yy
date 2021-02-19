%{
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <string>
#include <sstream>
#include <fstream>
#include <optional>
#include <functional>
#include <numeric>
#include "CTLLexer.ll.h"
#include "CTL.h"
#include "rgmedd5.h"
#include "parallel.h"

// #define PERFORMANCECTL 1
// #define PERFORMANCE 1

using namespace std;
using namespace ctlmdd;

// result of the formuls
static BaseFormula *g_parser_result;
// lexer instance
static FlexLexer* s_lexer;

void initialize_lexer(istringstream *p_iss) {
    assert(s_lexer == nullptr);
    s_lexer = new yyFlexLexer(p_iss);
}
void deinitialize_lexer() {
    assert(s_lexer != nullptr);
    delete s_lexer;
    s_lexer = nullptr;
}

extern const std::vector<size_t> *p_spot_ap_to_greatspn_ap_index;
extern const std::vector<Formula*> *p_greatspn_atomic_propositions;

void yyerror(const char *str) {
  cout<<"Parse error at \"" << s_lexer->YYText() << "\": " << str << "." << endl;
}
 
int yylex(void){
    int i = s_lexer->mmlex(); 
    return i;
}
extern int yyparse(void);
extern int output_flag;

// reverse the sign of an inequality (i.e. a < b -> b > a)
inline Inequality::op_type reverse_ineq_op(Inequality::op_type inop) {
    switch (inop) {
        case Inequality::IOP_MIN:       return Inequality::IOP_MAJ;
        case Inequality::IOP_MAJ:       return Inequality::IOP_MIN;
        case Inequality::IOP_MINEQ:     return Inequality::IOP_MAJEQ;
        case Inequality::IOP_MAJEQ:     return Inequality::IOP_MINEQ;
        default:                        return inop; // [not] equal
    }
}

inline Formula* fix_unquantified_ctlstar_formulas(Formula* f) {
    // if (f->isPathFormula()) { // typeid(*f) != typeid(QuantifiedFormula) && 
    //     /* Formula* negFormula = new LogicalFormula(f); */
    //     /* f = new QuantifiedFormula(f, QOP_ALWAYS); */
    //     f = ctlnew<QuantifiedFormula>(f, QOP_ALWAYS);
	// }
    return f;
}

//-----------------------------------------------------------------------------
%}
%union{
  float num;
  char *pVar;
  int place_id;
  int spot_id;
  int mpar_id;
  int transition_id;
  ctlmdd::PlaceTerm *term;
  ctlmdd::Formula *formula;
  ctlmdd::IntFormula *int_formula;
  ctlmdd::Inequality::op_type inop;
  std::vector<int>* place_id_list;
  std::vector<int>* transition_id_list;
}

%token <num> NUMBER
%token <pVar> VAR 
%token <mpar_id> MARK_PAR
%token <place_id> PLACE_ID
%token <transition_id> TRANSITION_ID
%token PROP_NAME SPOT_ACCEPT_ALL
%token PLUS MINUS TIMES DIV MINOR MAJOR MINOREQ MAJOREQ EQ NEQ 
%token OR XOR AND NOT IMPLY BIIMPLY POSSIBLY IMPOSSIBLY INVARIANT
%token HAS_DEADLOCK QUASI_LIVENESS STABLE_MARKING LIVENESS ONESAFE
%token LPARENT RPARENT TRUEv FALSEv LQPARENT RQPARENT 
%token DEADLOCK NDEADLOCK ENABLED BOUNDS COMMA
%token SIM DIF SHARP SEMICOLON
%right E A U ENABLED BOUNDS POSSIBLY IMPOSSIBLY INVARIANT // AX AF AG EX EF EG
%right X F G EX EF EG AX AF AG
%left IMPLY BIIMPLY
%left OR XOR PIPE
%left AND AMPERSAND
%right NOT
//%nonassoc SIM DIF
%nonassoc EQ MINOR MAJOR MINOREQ MAJOREQ NEQ
%left PLUS MINUS
%left TIMES DIV
%token LTLStart // switch to parse spot_expression

%type <int_formula> expression
%type <place_id_list> place_list
%type <transition_id_list> transition_list
%type <formula> atomic_prop
%type <inop> ineq_op
%type <formula> ctlstar_formula 
%type <formula> spot_expression

%start inizio

%%
inizio: expression opt_semicolon  { g_parser_result = $1; }
      // | ctl_formula opt_semicolon { parsed_comment = false; g_parser_result = $1; }
	  | LTLStart spot_expression  { g_parser_result = $2; } /* parse LTL atomic propositions generated by SPOT */
      | ctlstar_formula opt_semicolon { g_parser_result = fix_unquantified_ctlstar_formulas($1); }
      // | PROP_NAME VAR             { parsed_comment = true; g_parser_result = NULL; query_name = $2; free($2); }
      // | /*empty*/                 { parsed_comment = true; g_parser_result = NULL; }
      ;

opt_semicolon: /*nothing*/ | SEMICOLON ;


/** Boolean expression parser for SPOT's edge labels **/
spot_expression: spot_expression AMPERSAND spot_expression { $$ = ctlnew<LogicalFormula>($1, $3, LogicalFormula::CBF_AND);}
			   | spot_expression AND spot_expression       { $$ = ctlnew<LogicalFormula>($1, $3, LogicalFormula::CBF_AND);}
           	   | spot_expression PIPE spot_expression      { $$ = ctlnew<LogicalFormula>($1, $3, LogicalFormula::CBF_OR); }
               | spot_expression OR spot_expression        { $$ = ctlnew<LogicalFormula>($1, $3, LogicalFormula::CBF_OR); }
           	   | NOT spot_expression                       { $$ = ctlnew<LogicalFormula>($2); }
               | SPOT_ACCEPT_ALL                           { $$ = ctlnew<BoolLiteral>(true); }
               | NUMBER
               {
                    // Atomic proposition index must be present in the corresponding array
                    assert(p_greatspn_atomic_propositions != nullptr);
                    if ($1 > p_spot_ap_to_greatspn_ap_index->size() || $1 < 0) {
                        throw "ERROR: Atomic Proposition index is not valid."; // "
                    }
                    size_t ap_index = (*p_spot_ap_to_greatspn_ap_index)[$1];
                    $$ = (*p_greatspn_atomic_propositions)[ap_index];
                    $$->addOwner();
               }
               ;

/** CTLStar (CTL*) Expression **/
ctlstar_formula: atomic_prop                       { $$ = $1; }
             | LPARENT ctlstar_formula RPARENT     { $$ = $2; }
             | ctlstar_formula AND ctlstar_formula { $$ = ctlnew<LogicalFormula>($1,$3, LogicalFormula::CBF_AND); }
             | ctlstar_formula OR ctlstar_formula  { $$ = ctlnew<LogicalFormula>($1,$3, LogicalFormula::CBF_OR); }
             | NOT ctlstar_formula                 { $$ = ctlnew<LogicalFormula>($2); }
             | ctlstar_formula IMPLY ctlstar_formula { $$ = ctlnew<LogicalFormula>($1,$3, LogicalFormula::CBF_IMPLY); }
             | ctlstar_formula BIIMPLY ctlstar_formula { m_assert(false, "TODO: BIIMPLY"); }
             | POSSIBLY ctlstar_formula            { $$ = ctlnew<Reachability>($2, Reachability::RPT_POSSIBILITY); }
             | IMPOSSIBLY ctlstar_formula          { $$ = ctlnew<Reachability>($2, Reachability::RPT_IMPOSSIBILITY); }
             | INVARIANT ctlstar_formula           { $$ = ctlnew<Reachability>($2, Reachability::RPT_INVARIANTLY); }
             | A ctlstar_formula                   { $$ = ctlnew<QuantifiedFormula>($2, QOP_ALWAYS); }
             | E ctlstar_formula                   { $$ = ctlnew<QuantifiedFormula>($2, QOP_EXISTS); }
             | X ctlstar_formula                   { $$ = ctlnew<TemporalFormula>($2, POT_NEXT); }
             | G ctlstar_formula                   { $$ = ctlnew<TemporalFormula>($2, POT_GLOBALLY); }
             | F ctlstar_formula                   { $$ = ctlnew<TemporalFormula>($2, POT_FUTURE); }
             | ctlstar_formula U ctlstar_formula   { $$ = ctlnew<TemporalFormula>($1, $3); }
             | LQPARENT ctlstar_formula U ctlstar_formula RQPARENT { $$ = ctlnew<TemporalFormula>($2, $4); }
             /* syntactic sugar */
             | EX ctlstar_formula                  { $$ = ctlnew<QuantifiedFormula>(ctlnew<TemporalFormula>($2, POT_NEXT), QOP_EXISTS); }
             | EG ctlstar_formula                  { $$ = ctlnew<QuantifiedFormula>(ctlnew<TemporalFormula>($2, POT_GLOBALLY), QOP_EXISTS); }
             | EF ctlstar_formula                  { $$ = ctlnew<QuantifiedFormula>(ctlnew<TemporalFormula>($2, POT_FUTURE), QOP_EXISTS); }
             | AX ctlstar_formula                  { $$ = ctlnew<QuantifiedFormula>(ctlnew<TemporalFormula>($2, POT_NEXT), QOP_ALWAYS); }
             | AG ctlstar_formula                  { $$ = ctlnew<QuantifiedFormula>(ctlnew<TemporalFormula>($2, POT_GLOBALLY), QOP_ALWAYS); }
             | AF ctlstar_formula                  { $$ = ctlnew<QuantifiedFormula>(ctlnew<TemporalFormula>($2, POT_FUTURE), QOP_ALWAYS); }
             /* global properties */
             | HAS_DEADLOCK                        { $$ = ctlnew<GlobalProperty>(GPT_HAS_DEADLOCK); }
             | QUASI_LIVENESS                      { $$ = ctlnew<GlobalProperty>(GPT_QUASI_LIVENESS); }
             | STABLE_MARKING                      { $$ = ctlnew<GlobalProperty>(GPT_STABLE_MARKING); }
             | LIVENESS                            { $$ = ctlnew<GlobalProperty>(GPT_LIVENESS); }
             | ONESAFE                             { $$ = ctlnew<GlobalProperty>(GPT_ONESAFE); }
             ;

atomic_prop: NDEADLOCK                       { $$ = ctlnew<Deadlock>(false); }
           | DEADLOCK                        { $$ = ctlnew<Deadlock>(true); }
           | TRUEv                           { $$ = ctlnew<BoolLiteral>(true); }
           | FALSEv                          { $$ = ctlnew<BoolLiteral>(false); }
           | ENABLED LPARENT transition_list RPARENT  { $$ = ctlnew<Fireability>($3); delete $3; }
           | expression ineq_op expression   { $$ = make_inequality($1, $2, $3); }
           ;

ineq_op: EQ        { $$ = Inequality::IOP_EQ; }
       | MINOR     { $$ = Inequality::IOP_MIN; }
       | MINOREQ   { $$ = Inequality::IOP_MINEQ; }
       | MAJOR     { $$ = Inequality::IOP_MAJ; }
       | MAJOREQ   { $$ = Inequality::IOP_MAJEQ; }
       | NEQ       { $$ = Inequality::IOP_NEQ; }
       ;

place_list: opt_sharp PLACE_ID                  { $$ = new std::vector<int>(); $$->push_back($2); }
          | place_list COMMA opt_sharp PLACE_ID { $$ = $1; $$->push_back($4); }

transition_list: /* nothing */                       { $$ = new std::vector<int>(); }
               | TRANSITION_ID                       { $$ = new std::vector<int>(); $$->push_back($1); }
               | transition_list COMMA TRANSITION_ID { $$ = $1; $$->push_back($3); }

expression: LPARENT expression RPARENT        { $$ = $2;}
          | opt_sharp PLACE_ID                { $$ = ctlnew<PlaceTerm>(1, $2, PlaceTerm::EOP_TIMES); }
          | BOUNDS LPARENT place_list RPARENT { $$ = ctlnew<BoundOfPlaces>($3); delete $3; }
          | NUMBER                            { $$ = ctlnew<IntLiteral>($1); }
          | MARK_PAR                          { $$ = ctlnew<IntLiteral>(tabmp[$1].mark_val); }
          | MINUS expression  %prec NOT       { $$ = make_expression(ctlnew<IntLiteral>(0), IntFormula::EOP_MINUS, $2); }
          | expression TIMES expression       { $$ = make_expression($1, IntFormula::EOP_TIMES, $3); }
          | expression DIV expression         { $$ = make_expression($1, IntFormula::EOP_DIV, $3); }
          | expression PLUS expression        { $$ = make_expression($1, IntFormula::EOP_PLUS, $3); }
          | expression MINUS expression       { $$ = make_expression($1, IntFormula::EOP_MINUS, $3); }
          ;


opt_sharp : /*nothing*/ | SHARP;


%%

//-----------------------------------------------------------------------------

// Create an Inequality* object, with some optimizations for the special cases
AtomicProposition* make_inequality(IntFormula* e1, Inequality::op_type op, IntFormula* e2) {
    bool e1const = (typeid(*e1) == typeid(IntLiteral));
    bool e2const = (typeid(*e2) == typeid(IntLiteral));
    bool e1term = (typeid(*e1) == typeid(PlaceTerm));
    bool e2term = (typeid(*e2) == typeid(PlaceTerm));
    // constant <op> constant   ->   can be replaced with true/false
    if (e1const && e2const) {
        float val1 = ((IntLiteral*)e1)->getConstant();
        float val2 = ((IntLiteral*)e2)->getConstant();
        e1->removeOwner();
        e2->removeOwner();
        bool result;
        switch (op) {
            case Inequality::IOP_MIN:     result = val1 < val2;    break;
            case Inequality::IOP_MINEQ:   result = val1 <= val2;   break;
            case Inequality::IOP_MAJ:     result = val1 > val2;    break;
            case Inequality::IOP_MAJEQ:   result = val1 >= val2;   break;
            case Inequality::IOP_EQ:      result = val1 == val2;   break;
            case Inequality::IOP_NEQ:     result = val1 != val2;   break;
            case Inequality::IOP_SIM:     result = val1 == val2;   break;
            case Inequality::IOP_DIF:     result = val1 != val2;   break;
            default: throw;
        }
        return ctlnew<BoolLiteral>(result);
    }
    // constant <op> expression  ->  reverse the operator and build an inequality with constant
    else if (e1const) {
        float val1 = ((IntLiteral*)e1)->getConstant();
        e1->removeOwner();
        return ctlnew<Inequality>(reverse_ineq_op(op), e2, val1);
    }
    // expression <op> constant  ->  inequality with constant
    else if (e2const) {
        float val2 = ((IntLiteral*)e2)->getConstant();
        e2->removeOwner();
        return ctlnew<Inequality>(op, e1, val2);
    }
    // remaining case:  expression <op> expression
    // Use SIM and DIF if the two expressions are simple terms.    
    if (e1term && e2term) {
        if (op == Inequality::IOP_EQ)
            op = Inequality::IOP_SIM;
        else if (op == Inequality::IOP_NEQ)
            op = Inequality::IOP_DIF;
    }
    return ctlnew<Inequality>(op, e1, e2);

}

//-----------------------------------------------------------------------------

IntFormula* make_expression(IntFormula* e1, IntFormula::op_type op, IntFormula* e2) {
    bool e1const = (typeid(*e1) == typeid(IntLiteral));
    bool e2const = (typeid(*e2) == typeid(IntLiteral));
    bool e2term = (typeid(*e2) == typeid(PlaceTerm));
    // Terms are constants -> combine them directly
    if (e1const && e2const) {
        float result;
        float val1 = ((IntLiteral*)e1)->getConstant();
        float val2 = ((IntLiteral*)e2)->getConstant();
        switch (op) {
            case IntFormula::EOP_TIMES:   result = val1 * val2;     break;
            case IntFormula::EOP_DIV:     result = val1 / val2;     break;
            case IntFormula::EOP_PLUS:    result = val1 + val2;     break;
            case IntFormula::EOP_MINUS:   result = val1 - val2;     break;
            default: throw;
        }
        e1->removeOwner();
        e2->removeOwner();
        return ctlnew<IntLiteral>(result);
    }
    // <constant> <*/> <PlaceTerm>  ->  combine into a single PlaceTerm
    else if (e1const && e2term) {
        if (op == IntFormula::EOP_TIMES || op == IntFormula::EOP_DIV) {
            int variable = ((PlaceTerm*)e2)->getVariable();
            float coeff = ((PlaceTerm*)e2)->getCoeff();
            float val1 = ((IntLiteral*)e1)->getConstant();
            assert(coeff == 1);
            e1->removeOwner();
            e2->removeOwner();
            return ctlnew<PlaceTerm>(val1, variable, op);
        }
    }
    // Otherwise, create an IntFormula* object
    return ctlnew<IntExpression>(e1, e2, op);
}

//-----------------------------------------------------------------------------

BaseFormula* parse_formula() {
    assert(g_parser_result == nullptr);
    yyparse();
    BaseFormula* f = g_parser_result;
    g_parser_result = nullptr;
    return f;
}

//-----------------------------------------------------------------------------
