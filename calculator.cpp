#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define OK 1
#define ERROR 0

typedef int Status;

// 用来存储用户输入的向量或多项式或表达式的字符串
char a[1000];
char b[1000];

// 表示计算器开关状态的变量
int calculator_on = 1;

// 表示用户功能选择的变量
int selection;

// 用来清空输入缓冲区的辅助字符
char c;

char oprt[8] = {'+', '-', '*', '/', '^', '(', ')', '#'};

char priority[8][8] = {{'>', '>', '<', '<', '<', '<', '>', '>'},
                       {'>', '>', '<', '<', '<', '<', '>', '>'},
                       {'>', '>', '>', '>', '<', '<', '>', '>'},
                       {'>', '>', '>', '>', '<', '<', '>', '>'},
                       {'>', '>', '>', '>', '>', '<', '>', '>'},
                       {'<', '<', '<', '<', '<', '<', '=', '?'},
                       {'?', '?', '?', '?', '?', '?', '?', '?'},
                       {'<', '<', '<', '<', '<', '<', '?', '='}};

// *********************************************用来计算向量的顺序表Vector*********************************************
typedef struct {
    double *cmpt;
    int dim;
} Vector;

Status InitVector(Vector &V) {
    // 给向量分配存储空间，分配失败则报错
    V.cmpt = (double *)malloc(sizeof(double));
    if (!V.cmpt) {
        exit(OVERFLOW);
    }

    V.cmpt[0] = 0;
    V.dim = 0;

    return OK;
}

Status CreatVector(Vector &V, char *v) {
    // 给向量分配存储空间，分配失败则报错
    V.cmpt = (double *)malloc(V.dim * sizeof(double));
    if (!V.cmpt){
        exit(OVERFLOW);
    }

    int i = 0;
    char *space_addr = v;

    // 每检测到一个空格时，读取空格后面的数字
    if (v[0] == '[') {
        sscanf(v + 1, "%lf", &V.cmpt[0]);
    } else {
        sscanf(v, "%lf", &V.cmpt[0]);
    }
    while (strchr(space_addr + 1, ' ')) {
        space_addr = strchr(space_addr + 1, ' ');
        sscanf(space_addr + 1, "%lf", &V.cmpt[++i]);
    }

    return OK;
}

Status DestroyVector(Vector &V) {
    // 释放空间
    free(V.cmpt);
    V.cmpt = nullptr;

    Vector *p = &V;
    p = nullptr;

    return OK;
}

Vector CalcVector(Vector A, Vector B, char opr) {
    Vector result;

    InitVector(result);

    // 给向量分配存储空间，分配失败则报错
    result.cmpt = (double *)realloc(result.cmpt, A.dim * sizeof(double));
    if (!result.cmpt) {
        exit(OVERFLOW);
    }

    result.dim = A.dim;

    // 根据运算符，使向量result等于向量A与B的和或差
    for (int i = 0; i < A.dim; ++i) {
        switch (opr) {
            case '+':
                result.cmpt[i] = A.cmpt[i] + B.cmpt[i];
                break;
            case '-':
                result.cmpt[i] = A.cmpt[i] - B.cmpt[i];
                break;
            default:;
        }
    }

    return result;
}

double CalcAngleBetweenVector(Vector A, Vector B) {
    double a_norm = 0, b_norm = 0;
    double dot_product = 0;

    for (int i = 0; i < A.dim; ++i) {
        a_norm += A.cmpt[i] * A.cmpt[i];
        b_norm += B.cmpt[i] * B.cmpt[i];
        dot_product += A.cmpt[i] * B.cmpt[i];
    }

    a_norm = sqrt(a_norm);
    b_norm = sqrt(b_norm);

    // 将上述结果代入向量夹角的计算公式
    return dot_product / (a_norm * b_norm);
}

void PrintVector(Vector V) {
    putchar('[');

    for (int i = 0; i < V.dim - 1; ++i) {
        printf("%.2f ", V.cmpt[i]);
    }

    printf("%.2f]", V.cmpt[V.dim - 1]);
}

// ***************************************用来计算一元多项式的顺序表SeqPolynomial***************************************
typedef struct {
    double coef;
    int expn;
    int num_of_terms;
} SeqTerm, *SeqPolynomial;

Status InitSeqPoly(SeqPolynomial &P) {
    // 给多项式分配存储空间，分配失败则报错
    P = (SeqPolynomial)malloc(sizeof(SeqTerm));
    if (!P) {
        exit(OVERFLOW);
    }

    P[0].coef = 0;
    P[0].expn = 0;
    P->num_of_terms = 0;

    return OK;
}

Status DestroySeqPoly(SeqPolynomial &P) {
    // 释放空间
    free(P);
    P = nullptr;

    return OK;
}

int LocateSeqTerm(SeqPolynomial P, SeqTerm t, int *loc) {
    // 遍历顺序表，找到第一个次数不超过待插项的项
    for (int i = 0; i < P->num_of_terms; ++i) {
        if (P[i].expn == t.expn) {
            *loc = i;
            return 1;
        }
        if (P[i].expn < t.expn) {
            *loc = i;
            return 0;
        }
    }

    *loc = P->num_of_terms;
    return 0;
}

Status DeleteSeqTerm(SeqPolynomial &P, int loc) {
    if (P->num_of_terms == 1) {
        P->num_of_terms = 0;
        return OK;
    } else if (!loc) {
        P[1].num_of_terms = P->num_of_terms;
    }

    for (int i = loc; i < P->num_of_terms - 1; ++i) {
        P[i] = P[i + 1];
    }

    P = (SeqPolynomial)realloc(P, --P->num_of_terms * sizeof(SeqTerm));

    return OK;
}

Status InsertSeqTerm(SeqPolynomial &P, SeqTerm t, char opr) {
    if (!P->num_of_terms) {
        P[0] = t;
        P->num_of_terms = 1;
        return OK;
    }

    int loc = 0;

    if (LocateSeqTerm(P, t, &loc)) {
        switch(opr) {
            case '+':
                P[loc].coef += t.coef;
                break;
            case '-':
                P[loc].coef -= t.coef;
                break;
            default:;
        }

        if (!P[loc].coef) {
            DeleteSeqTerm(P, loc);
        }
    } else {
        P = (SeqPolynomial)realloc(P, ++P->num_of_terms * sizeof(SeqTerm));
        if (!P) {
            exit(OVERFLOW);
        }

        for (int i = P->num_of_terms - 1; i > loc; --i) {
            P[i] = P[i - 1];
        }

        P[loc] = t;
        if (opr == '-') {
            P[loc].coef *= -1;
        }

        if (!loc) {
            P->num_of_terms = P[1].num_of_terms;
        }
    }

    return OK;
}

Status CreatSeqPoly(SeqPolynomial &P, char *p) {
    // 给多项式分配存储空间，分配失败则报错
    int i = 0;
    char *space_addr = p;

    SeqTerm t;

    // 每检测到一个空格时，读取空格后面的数字
    sscanf(p, "%lf", &t.coef);
    while (strchr(space_addr + 1, ' ')) {
        space_addr = strchr(space_addr + 1, ' ');
        if (++i % 2) {
            sscanf(space_addr + 1, "%d", &t.expn);
            InsertSeqTerm(P, t, '+');
        } else {
            sscanf(space_addr + 1, "%lf", &t.coef);
        }
    }

    return OK;
}

SeqPolynomial CalcSeqPoly(SeqPolynomial A, SeqPolynomial B, char opr) {
    SeqPolynomial result = nullptr;

    InitSeqPoly(result);

    if (opr == '+' || opr == '-') {
        for (int i = 0; i < A->num_of_terms; ++i) {
            InsertSeqTerm(result, A[i], '+');
        }

        for (int i = 0; i < B->num_of_terms; ++i) {
            InsertSeqTerm(result, B[i], opr);
        }
    } else if (opr == '*') {
        SeqTerm t;

        for (int i = 0; i < A->num_of_terms; ++i) {
            for (int j = 0; j < B->num_of_terms; ++j) {
                t.coef = A[i].coef * B[j].coef;
                t.expn = A[i].expn + B[j].expn;
                InsertSeqTerm(result, t, '+');
            }
        }
    }

    return result;
}

SeqPolynomial CalcDerivSeqPoly(SeqPolynomial A, int deriv_order) {
    SeqPolynomial result = nullptr;

    InitSeqPoly(result);

    result = A;

    if (!deriv_order) {
        return result;
    }

    for (int i = 0; i < result->num_of_terms; ++i) {
        result[i].coef *= result[i].expn;
        if (!result[i].coef) {
            if (result->num_of_terms == 1) {
                result->num_of_terms = 0;
                return result;
            } else {
                DeleteSeqTerm(result, i--);
                continue;
            }
        }
        --result[i].expn;
    }

    return CalcDerivSeqPoly(result, deriv_order - 1);
}

void PrintSeqPoly(SeqPolynomial P) {
    if (!P->num_of_terms) {
        printf("0");
        return;
    }

    // 输出第一项
    if (!P[0].expn) {
        printf("%.2f", P[0].coef);
    } else if (P[0].expn == 1 && P[0].coef == 1) {
        printf("x");
    } else if (P[0].expn == 1 && P[0].coef == -1) {
        printf("-x");
    } else if (P[0].expn == 1) {
        printf("%.2fx", P[0].coef);
    } else if (P[0].coef == 1) {
        printf("x ^ %d", P[0].expn);
    } else if (P[0].coef == -1) {
        printf("-x ^ %d", P[0].expn);
    } else {
        printf("%.2fx ^ %d", P[0].coef, P[0].expn);
    }

    // 输出剩余项
    for (int i = 1; i < P->num_of_terms; ++i) {
        if (!P[i].expn) {
            if (P[i].coef < 0) {
                printf(" - %.2f",  -P[i].coef);
            } else {
                printf(" + %.2f", P[i].coef);
            }
        } else if (P[i].expn == 1 && P[i].coef == 1) {
            printf(" + x");
        } else if (P[i].expn == 1 && P[i].coef == -1) {
            printf(" - x");
        } else if (P[i].expn == 1) {
            if (P[i].coef < 0) {
                printf(" - %.2fx", -P[i].coef);
            } else {
                printf(" + %.2fx", P[i].coef);
            }
        } else if (P[i].coef == 1) {
            printf(" + x ^ %d", P[i].expn);
        } else if (P[i].coef == -1) {
            printf(" - x ^ %d", P[i].expn);
        } else {
            if (P[i].coef < 0) {
                printf(" - %.2fx ^ %d", -P[i].coef, P[i].expn);
            } else {
                printf(" + %.2fx ^ %d", P[i].coef, P[i].expn);
            }
        }
    }
}

// ***************************************用来计算一元多项式的链表LinkPolynomial***************************************
typedef struct LinkTerm {
    double coef;
    int expn;
    struct LinkTerm *next;
} LinkTerm, *LinkPolynomial;

Status InitLinkPoly(LinkPolynomial &P) {
    // 给多项式分配存储空间，分配失败则报错
    P = (LinkPolynomial)malloc(sizeof(LinkTerm));
    if (!P) {
        exit(OVERFLOW);
    }

    P->next = nullptr;

    return OK;
}

Status DestroyLinkPoly(LinkPolynomial &P) {
    LinkTerm *pre = P;

    // 遍历各结点，分别释放空间
    while (P) {
        P = P->next;
        free(pre);
        pre = P;
    }

    return OK;
}

int LocateLinkTerm(LinkPolynomial P, LinkTerm t, LinkTerm *&pre) {
    pre = P;
    LinkTerm *p;
    p = P->next;

    // 遍历链表，找到第一个次数不超过待插项的项的前驱
    while (p && p->expn > t.expn) {
        pre = p;
        p = p->next;
    }

    if (p && p->expn == t.expn) {
        return 1;
    } else {
        return 0;
    }
}

Status InsertLinkTerm(LinkPolynomial &P, LinkTerm t, char opr) {
    LinkTerm *pre = nullptr;

    if (LocateLinkTerm(P, t, pre)) {
        switch(opr) {
            case '+':
                pre->next->coef += t.coef;
                break;
            case '-':
                pre->next->coef -= t.coef;
                break;
            default:;
        }

        if (!pre->next->coef) {
            LinkTerm *p = pre->next;
            pre->next = p->next;
            free(p);
        }
    } else {
        LinkTerm *T = (LinkTerm *)malloc(sizeof(LinkTerm));
        *T = t;
        if (opr == '-') {
            T->coef *= -1;
        }

        T->next = pre->next;
        pre->next = T;
    }

    return OK;
}

Status CreatLinkPoly(LinkPolynomial &P, char *p) {
    char *space_addr = p;

    LinkTerm t;

    while (space_addr) {
        if (space_addr == p) {
            sscanf(space_addr, "%lf", &t.coef);
            space_addr = strchr(space_addr, ' ');
        } else {
            sscanf(space_addr + 1, "%lf", &t.coef);
            space_addr = strchr(space_addr + 1, ' ');
        }

        sscanf(space_addr + 1, "%d", &t.expn);
        space_addr = strchr(space_addr + 1, ' ');

        InsertLinkTerm(P, t, '+');
    }

    return OK;
}

LinkPolynomial CalcLinkPoly(LinkPolynomial A, LinkPolynomial B, char opr) {
    LinkPolynomial result = nullptr;

    InitLinkPoly(result);

    LinkTerm *t = (LinkTerm *)malloc(sizeof(LinkTerm));

    if (opr == '+' || opr == '-') {
        t = A->next;

        while (t) {
            InsertLinkTerm(result, *t, '+');
            t = t->next;
        }

        t = B->next;

        while (t) {
            InsertLinkTerm(result, *t, opr);
            t = t->next;
        }
    } else if (opr == '*') {
        LinkTerm *pa = A->next;
        LinkTerm *pb = B->next;

        while (pa) {
            while (pb) {
                t->coef = pa->coef * pb->coef;
                t->expn = pa->expn + pb->expn;
                InsertLinkTerm(result, *t, '+');
                pb = pb->next;
            }
            pa = pa->next;
            pb = B->next;
        }
    }

    return result;
}

LinkPolynomial CalcDerivLinkPoly(LinkPolynomial A, int deriv_order) {
    LinkPolynomial result = nullptr;

    InitLinkPoly(result);

    result = A;

    if (!deriv_order) {
        return result;
    }

    LinkTerm *p = result;

    while (p->next) {
        p->next->coef *= p->next->expn;
        if (!p->next->coef) {
            LinkTerm *ptr = p->next;
            p->next = ptr->next;
            free(ptr);
            continue;
        }

        --p->next->expn;
        p = p->next;
    }

    return CalcDerivLinkPoly(result, deriv_order - 1);
}

void PrintLinkPoly(LinkPolynomial P) {
    if (!P->next) {
        printf("0");
        return;
    }

    LinkTerm *p = P->next;

    // 输出第一项
    if (!p->expn) {
        printf("%.2f", p->coef);
    } else if (p->expn == 1 && p->coef == 1) {
        printf("x");
    } else if (p->expn == 1 && p->coef == -1) {
        printf("-x");
    } else if (p->expn == 1) {
        printf("%.2fx", p->coef);
    } else if (p->coef == 1) {
        printf("x ^ %d", p->expn);
    } else if (p->coef == -1) {
        printf("-x ^ %d", p->expn);
    } else {
        printf("%.2fx ^ %d", p->coef, p->expn);
    }

    p = p->next;

    // 输出剩余项
    while (p) {
        if (!p->expn) {
            if (p->coef < 0) {
                printf(" - %.2f", -p->coef);
            } else {
                printf(" + %.2f", p->coef);
            }
        } else if (p->expn == 1 && p->coef == 1) {
            printf(" + x");
        } else if (p->expn == 1 && p->coef == -1) {
            printf(" - x");
        } else if (p->expn == 1) {
            if (p->coef < 0) {
                printf(" - %.2fx", -p->coef);
            } else {
                printf(" + %.2fx", p->coef);
            }
        } else if (p->coef == 1) {
            printf(" + x ^ %d", p->expn);
        } else if (p->coef == -1) {
            printf(" - x ^ %d", p->expn);
        } else {
            if (p->coef < 0) {
                printf(" - %.2fx ^ %d", -p->coef, p->expn);
            } else {
                printf(" + %.2fx ^ %d", p->coef, p->expn);
            }
        }

        p = p->next;
    }
}

// ******************************************用来进行表达式求值的顺序栈Stack******************************************
#define STACK_INIT_SIZE 10 // 存储空间初始分配量
#define STACK_INCREMENT 10 // 存储空间分配增量

typedef struct {
    double *base;
    double *top;
    int stacksize;
} Stack_OPND;

typedef struct {
    char *base;
    char *top;
    int stacksize;
} Stack_OPTR;

Status InitStack_OPND(Stack_OPND &S) {
    // 构造一个空栈
    S.base = (double *)malloc(STACK_INIT_SIZE * sizeof(double));
    if (!S.base) {
        exit(OVERFLOW);
    }

    S.top = S.base;
    S.stacksize = STACK_INIT_SIZE;

    return OK;
}

Status InitStack_OPTR(Stack_OPTR &S) {
    // 构造一个空栈
    S.base = (char *)malloc(STACK_INIT_SIZE * sizeof(char));
    if (!S.base) {
        exit(OVERFLOW);
    }

    S.top = S.base;
    S.stacksize = STACK_INIT_SIZE;

    return OK;
}

Status Push_OPND(Stack_OPND &S, double e) {
    // 若栈满，则追加存储空间
    if (S.top - S.base >= S.stacksize) {
        S.base = (double *)realloc(S.base, (S.stacksize + STACK_INCREMENT) * sizeof(double));
        if (!S.base) {
            exit(OVERFLOW);
        }
        S.top =  S.base + S.stacksize;
        S.stacksize += STACK_INCREMENT;
    }

    // 插入e为新的栈顶元素
    *S.top++ = e;

    return OK;
}

Status Push_OPTR(Stack_OPTR &S, char e) {
    // 若栈满，则追加存储空间
    if (S.top - S.base >= S.stacksize) {
        S.base = (char *)realloc(S.base, (S.stacksize + STACK_INCREMENT) * sizeof(char));
        if (!S.base) {
            exit(OVERFLOW);
        }
        S.top =  S.base + S.stacksize;
        S.stacksize += STACK_INCREMENT;
    }

    // 插入e为新的栈顶元素
    *S.top++ = e;

    return OK;
}

Status Pop_OPND(Stack_OPND &S, double &e) {
    // 若栈不为空，删除栈顶元素，并用e返回其值
    if (S.top == S.base) {
        return ERROR;
    }

    e = *--S.top;

    return OK;
}

Status Pop_OPTR(Stack_OPTR &S, char &e) {
    // 若栈不为空，删除栈顶元素，并用e返回其值
    if (S.top == S.base) {
        return ERROR;
    }

    e = *--S.top;

    return OK;
}

double GetTop_OPND(Stack_OPND S) {
    // 若栈不为空，返回栈顶元素，否则返回0
    if (S.top == S.base) {
        return ERROR;
    }

    return *(S.top - 1);
}

char GetTop_OPTR(Stack_OPTR S) {
    // 若栈不为空，返回栈顶元素，否则返回0
    if (S.top == S.base) {
        return ERROR;
    }

    return *(S.top - 1);
}

bool ExpressionIsInvalid(char *str) {
    int len = strlen(str);

    for (int i = 0; i < len; ++i) {
        if (!strchr(oprt, str[i]) && str[i] != ' ' && str[i] != '.' && str[i] != '_' && (str[i] < 48 || str[i] > 57) && (str[i] < 65 || str[i] > 90) && (str[i] < 97 || str[i] > 122)) {
            return true;
        }
    }

    return false;
}

char * VariableContained(char *str) {
    int len = strlen(str);

    for (int i = 0; i < len; ++i) {
        if ((str[i] >= 65 && str[i] <= 90) || (str[i] >= 97 && str[i] <= 122) || str[i] == '_') {
            int n = 1;
            while (i + n < len && !strchr(oprt, str[i + n])) {
                ++n;
            }

            char var[n];
            strncpy(var, str + i, n);

            return var;
        }
    }

    return nullptr;
}

void ReplaceVariable(char *str, char *var, char *value) {
    int len_var = strlen(var);
    int len_value = strlen(value);

    char *p = str;
    while (p < str + strlen(str) && strstr(p, var)) {
        p = strstr(p, var);
        if (len_var >= len_value) {
            for (int i = 0; i < len_value; ++i) {
                p[i] = value[i];
            }
            for (int i = 0; i < (str + strlen(str)) - (p + len_var); ++i) {
                p[i + len_value] = p[i + len_var];
            }
        } else {
            for (int i = str + strlen(str) - p + len_value - len_var - 1; i >= len_value; --i) {
                p[i] = p[i - (len_value - len_var)];
            }
            for (int i = 0; i < len_value; ++i) {
                p[i] = value[i];
            }
        }

        str[strlen(str) + len_value - len_var] = '\0';

        p += len_value;
    }
}

char Precede(char a, char b) {
    int a_index, b_index;

    for (int i = 0; i < 8; ++i) {
        if (a == oprt[i]) {
            a_index = i;
        }
        if (b == oprt[i]) {
            b_index = i;
        }
    }

    return priority[a_index][b_index];
}

double Operate(double a, char theta, double b) {
    switch (theta) {
        case '+':
            return (a + b);
            break;
        case '-':
            return (a - b);
            break;
        case '*':
            return (a * b);
            break;
        case '/':
            return (a / b);
            break;
        case '^':
            return pow(a, b);
            break;
        default:;
    }
}

double CalcExpression(char *str) {
    int len = strlen(str);
    str[len] = '#';
    str[++len] = '\0';

    Stack_OPND OPND;
    Stack_OPTR OPTR;

    InitStack_OPND(OPND);
    InitStack_OPTR(OPTR);

    Push_OPTR(OPTR, '#');

    int NegetiveSign = 1;
    bool DetectDecimal = false;
    int count = 0;

    for (int i = 0; i < len; ++i) {
        if (!strchr(oprt, str[i])) {
            // 如果读到的字符是数字或小数点
            if (i == 0) {
                Push_OPND(OPND, str[i] - '0');
            } else if (str[i] == '.') {
                DetectDecimal = true;
            } else if (!strchr(oprt, str[i - 1])) {
                double x;

                if (DetectDecimal) {
                    // 若正在读入的是小数点后的数字
                    ++count;
                    Pop_OPND(OPND, x);
                    Push_OPND(OPND, x + NegetiveSign * pow(10, -count) * (str[i] - '0'));
                } else {
                    // 若正在读入的是小数点前的数字
                    Pop_OPND(OPND, x);
                    Push_OPND(OPND, 10 * x + NegetiveSign * (str[i] - '0'));
                }
            } else {
                Push_OPND(OPND, NegetiveSign * (str[i] - '0'));
            }
        } else {
            // 如果读到的字符是运算符
            DetectDecimal = false;
            count = 0;

            if (str[i] == '-' && (i == 0 || (str[i - 1] != ')' && strchr(oprt, str[i - 1])))) {
                // 若读取到的'-'为负号
                if (str[i - 1] == '-' && (i == 1 || (str[i - 2] != ')' && strchr(oprt, str[i - 2])))) {
                    NegetiveSign *= -1;
                } else {
                    NegetiveSign = -1;
                }
            } else if (str[i] == '-' && (str[i - 1] == '-' && (i == 1 || (str[i - 2] != ')' && strchr(oprt, str[i - 2]))))) {
                NegetiveSign = 1;
            } else {
                NegetiveSign = 1;

                switch (Precede(GetTop_OPTR(OPTR), str[i])) {
                    case '<':
                        Push_OPTR(OPTR, str[i]);
                        break;
                    case '=':
                        char _c;
                        Pop_OPTR(OPTR, _c);
                        break;
                    case '>':
                        double _a, _b;
                        char theta;
                        Pop_OPTR(OPTR, theta);
                        Pop_OPND(OPND, _b);
                        Pop_OPND(OPND, _a);
                        Push_OPND(OPND, Operate(_a, theta, _b));
                        --i;
                        break;
                    default:;
                }
            }
        }
    }

    return GetTop_OPND(OPND);
}

// *************************************************函数储存的数据结构*************************************************
typedef struct {
    char name[20] = "NAME";
    char variable[10] = "VARIABLE";
    char body[100] = "BODY";
    int num_of_func = 0;
} Function;

Function function[50];

Status AddFunction(char *f) {
    ++function->num_of_func;

    char *p1 = strchr(f, '(');
    char *p2 = strchr(f, ')');

    strncpy(function[function->num_of_func - 1].name, f, p1 - f);
    function[function->num_of_func - 1].name[p1 - f] = '\0';
    strncpy(function[function->num_of_func - 1].variable, p1 + 1, p2 - p1 - 1);
    function[function->num_of_func - 1].variable[p2 - p1 - 1] = '\0';
    strcpy(function[function->num_of_func - 1].body, p2 + 2);
    function[function->num_of_func - 1].body[strlen(f) - (p2 - f) - 2] = '\0';

    return OK;
}

// *************************************************辅助功能的实现代码*************************************************
void Menu() {
    // 用户挑选要使用的功能
    printf("\nSelect a function:\n");
    printf("Vector calculation......0    Polynomial calculation(using SeqList)......1\n"
           "Polynomial calculation(using LinkList)......2    Evaluation of expression......3\n"
           "Function calculation......4\n");
    printf(">>> ");

    scanf("%d", &selection);
    while((c = (char)getchar()) != '\n' && c != EOF);

    // 用户输入不合规数字时，要求重新输入
    while (selection < 0 || selection > 4) {
        printf("\nThis is not a valid selection. Please select again.\n");
        printf("Vector calculation......0    Polynomial calculation(using SeqList)......1\n"
               "Polynomial calculation(using LinkList)......2    Evaluation of expression......3\n"
               "Function calculation......4\n");
        printf(">>> ");
        scanf("%d", &selection);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }
}

void CalcComplete() {
    // 计算完成后向用户询问是否继续使用计算器
    char c;

    printf("\nCalculation completed. Continue using calculator?\n");
    printf("Continue......0    Exit......1\n");
    printf(">>> ");

    scanf("%d", &calculator_on);
    while((c = (char)getchar()) != '\n' && c != EOF);

    while (calculator_on != 0 && calculator_on != 1) {
        printf("\nThis is not a valid selection. Please select again.\n");
        printf("Continue......0    Exit......1\n");
        printf(">>> ");
        scanf("%d", &calculator_on);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }

    if (calculator_on) {
        calculator_on = 0;
    } else {
        calculator_on = 1;
    }
}

int CountSpace(char *str) {
    // 数字符串中的空格数量
    int count = 0;
    unsigned long len = strlen(str);

    for (int i = 0; i < len; ++i) {
        if (str[i] == ' ')
            ++count;
    }

    return count;
}

void ClearSpace(char *str) {
    // 清除用户输入的表达式中的空格
    int count = 0;
    int len = strlen(str);

    for (int i = 0; i < len; ++i) {
        if (str[i] == ' ') {
            ++count;
        } else {
            str[i - count] = str[i];
        }
    }

    str[len - count] = '\0';
}

// ***********************************************计算器各功能的运行代码***********************************************
void VectorCalculation() {
    printf("\nPlease input a vector as A:\n");
    printf(">>> ");
    scanf("%[^\n]", a);
    while((c = (char)getchar()) != '\n' && c != EOF);

    printf("\nPlease input a vector as B:\n");
    printf(">>> ");
    scanf("%[^\n]", b);
    while((c = (char)getchar()) != '\n' && c != EOF);

    Vector A, B;

    A.dim = CountSpace(a) + 1;
    B.dim = CountSpace(b) + 1;

    while (A.dim != B.dim) {
        printf("\nThe dimension of A and B are not equal. Please input another vector as B:\n");
        printf(">>> ");
        scanf("%[^\n]", b);
        while((c = (char)getchar()) != '\n' && c != EOF);

        B.dim = CountSpace(b) + 1;
    }

    CreatVector(A, a);
    CreatVector(B, b);

    printf("\nCalculating...\n\n");

    // 输出计算结果
    printf("A = ");
    PrintVector(A);
    printf("\nB = ");
    PrintVector(B);
    printf("\nA + B = ");
    PrintVector(CalcVector(A, B, '+'));
    printf("\nA - B = ");
    PrintVector(CalcVector(A, B, '-'));
    printf("\ncos Θ = %.4f\n", CalcAngleBetweenVector(A, B));

    DestroyVector(A);
    DestroyVector(B);
}

void SeqPolyCalculation() {
    printf("\nPlease input a polynomial as A:\n");
    printf(">>> ");
    scanf("%[^\n]", a);
    while((c = (char)getchar()) != '\n' && c != EOF);

    while (!(CountSpace(a) % 2)) {
        printf("\nInvalid input. Please input another polynomial as A:\n");
        printf(">>> ");
        scanf("%[^\n]", a);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }

    printf("\nPlease input a polynomial as B:\n");
    printf(">>> ");
    scanf("%[^\n]", b);
    while((c = (char)getchar()) != '\n' && c != EOF);

    while (!(CountSpace(b) % 2)) {
        printf("\nInvalid input. Please input another polynomial as B:\n");
        printf(">>> ");
        scanf("%[^\n]", b);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }

    int deriv_order;

    printf("\nPlease input the derivation order of A:\n");
    printf(">>> ");
    scanf("%d", &deriv_order);
    while((c = (char)getchar()) != '\n' && c != EOF);

    while (deriv_order < 0) {
        printf("\nDerivation number should be non-negative. Please input another number:\n");
        printf(">>> ");
        scanf("%d", &deriv_order);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }

    SeqPolynomial A, B;

    InitSeqPoly(A);
    InitSeqPoly(B);

    CreatSeqPoly(A, a);
    CreatSeqPoly(B, b);

    printf("\nCalculating...\n\n");

    // 输出计算结果
    printf("A = ");
    PrintSeqPoly(A);
    printf("\nB = ");
    PrintSeqPoly(B);
    printf("\nA + B = ");
    PrintSeqPoly(CalcSeqPoly(A, B, '+'));
    printf("\nA - B = ");
    PrintSeqPoly(CalcSeqPoly(A, B, '-'));
    printf("\nA * B = ");
    PrintSeqPoly(CalcSeqPoly(A, B, '*'));
    printf("\nThe %dth order derivative of A = ", deriv_order);
    PrintSeqPoly(CalcDerivSeqPoly(A, deriv_order));
    putchar('\n');

    DestroySeqPoly(A);
    DestroySeqPoly(B);
}

void LinkPolyCalculation() {
    printf("\nPlease input a polynomial as A:\n");
    printf(">>> ");
    scanf("%[^\n]", a);
    while((c = getchar()) != '\n' && c != EOF);

    while (!(CountSpace(a) % 2)) {
        printf("\nInvalid input. Please input another polynomial as A:\n");
        printf(">>> ");
        scanf("%[^\n]", a);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }

    printf("\nPlease input a polynomial as B:\n");
    printf(">>> ");
    scanf("%[^\n]", b);
    while((c = (char)getchar()) != '\n' && c != EOF);

    while (!(CountSpace(b) % 2)) {
        printf("\nInvalid input. Please input another polynomial as B:\n");
        printf(">>> ");
        scanf("%[^\n]", b);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }

    int deriv_order;

    printf("\nPlease input the derivation order of A:\n");
    printf(">>> ");
    scanf("%d", &deriv_order);
    while((c = (char)getchar()) != '\n' && c != EOF);

    while (deriv_order < 0) {
        printf("\nDerivation number should be non-negative. Please input another number:\n");
        printf(">>> ");
        scanf("%d", &deriv_order);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }

    LinkPolynomial A, B;

    InitLinkPoly(A);
    InitLinkPoly(B);

    CreatLinkPoly(A, a);
    CreatLinkPoly(B, b);

    printf("\nCalculating...\n\n");

    // 输出计算结果
    printf("A = ");
    PrintLinkPoly(A);
    printf("\nB = ");
    PrintLinkPoly(B);
    printf("\nA + B = ");
    PrintLinkPoly(CalcLinkPoly(A, B, '+'));
    printf("\nA - B = ");
    PrintLinkPoly(CalcLinkPoly(A, B, '-'));
    printf("\nA * B = ");
    PrintLinkPoly(CalcLinkPoly(A, B, '*'));
    printf("\nThe %dth order derivative of A = ", deriv_order);
    PrintLinkPoly(CalcDerivLinkPoly(A, deriv_order));
    putchar('\n');

    DestroyLinkPoly(A);
    DestroyLinkPoly(B);
}

void EvaluateExpression() {
    printf("\nPlease input an expression:\n");
    printf(">>> ");
    scanf("%[^\n]", a);
    while((c = (char)getchar()) != '\n' && c != EOF);

    while (ExpressionIsInvalid(a)) {
        printf("\nInvalid input. Please input another expression:\n");
        printf(">>> ");
        scanf("%[^\n]", a);
        while((c = (char)getchar()) != '\n' && c != EOF);
    }

    ClearSpace(a);

    char *var = VariableContained(a);

    if (var) {
        char value[10];

        printf("\nPlease input the value of %s:\n", var);
        printf(">>> ");
        scanf("%s", value);
        while((c = (char)getchar()) != '\n' && c != EOF);

        ReplaceVariable(a, var, value);
    }

    printf("\nCalculating...\n\n");
    printf("Answer = %g\n", CalcExpression(a));
}

void FunctionCalculation() {
    printf("\nStart to input command lines:\n");
    while (a[0] != 'E') {
        scanf("%[^\n]", a);
        while((c = getchar()) != '\n' && c != EOF);

        if (a[0] == 'D') {
            char *f = a + 4;
            ClearSpace(f);
            AddFunction(f);
        } else if (a[0] == 'R') {
            char *f = a + 4;

            char func_name[20] = "NAME";
            char func_value[10] = "VALUE";
            char func_body[100] = "BODY";

            char *p1 = strchr(f, '(');
            char *p2 = strchr(f, ')');

            strncpy(func_name, f, p1 - f);
            func_name[p1 - f] = '\0';
            strncpy(func_value, p1 + 1, p2 - p1 - 1);
            func_value[p2 - p1 - 1] = '\0';

            int i = 0;
            while (strcmp(func_name, function[i].name) && i < function->num_of_func) {
                ++i;
            }

            strcpy(func_body, function[i].body);
            func_body[strlen(function[i].body)] = '\0';

            ReplaceVariable(func_body, function[i].variable, func_value);

            printf("%g\n", CalcExpression(func_body));
        }
    }
}

// ******************************************************主函数******************************************************
int main() {
    printf("Welcome to Easy Calculator!\n");

    while (calculator_on) {
        Menu();

        switch (selection) {
            case 0:
                VectorCalculation();
                break;
            case 1:
                SeqPolyCalculation();
                break;
            case 2:
                LinkPolyCalculation();
                break;
            case 3:
                EvaluateExpression();
                break;
            case 4:
                FunctionCalculation();
                break;
            default:;
        }

        CalcComplete();
    }

    return 0;
}