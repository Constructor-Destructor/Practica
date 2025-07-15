#include <iostream>
#include <vector>
#include <iomanip>
#include <locale>
#include <clocale>
#include <limits>
#include <cmath>
#include <sstream>

using namespace std;

const int MAX_YZEL = 1000;
const double MIN_DIFF = 1e-12;

// Вычисление таблицы разделённых разностей
vector<vector<double>> tabl_razdelen_razn(const vector<double>& x, const vector<double>& y)
{
    int n = x.size();
    vector<vector<double>> table(n, vector<double>(n));
    for (int i = 0; i < n; i++)
    {
        table[i][0] = y[i];
    }

    for (int j = 1; j < n; j++)
    {
        for (int i = 0; i < n - j; i++)
        {
            double denom = x[i + j] - x[i];
            if (abs(denom) < MIN_DIFF)
            {
                cerr << "Ошибка: x[" << i + j << "] и x[" << i << "] слишком близки. Деление на почти ноль.\n";
                exit(1);
            }
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / denom;
        }
    }

    return table;
}

// Вывод формулы полинома Ньютона
void printPolinom(const vector<double>& x, const vector<vector<double>>& table)
{
    int n = x.size();
    cout << "\n Пошаговое построение полинома Ньютона:\n";
    cout << fixed << setprecision(6);

    cout << "\n Расчёт разделённых разностей:\n";
    for (int i = 0; i < n; i++)
    {
        if (i == 0)
        {
            cout << "f[x0] = y[0] = " << table[0][0] << endl;
        }
        else {
            cout << "f[";
            for (int j = 0; j <= i; j++)
            {
                if (j > 0) cout << ",";
                cout << "x" << j;
            }
            cout << "] = (f[";
            for (int j = 1; j <= i; j++)
            {
                if (j > 1) cout << ",";
                cout << "x" << j;
            }
            cout << "] - f[";
            for (int j = 0; j < i; j++)
            {
                if (j > 0) cout << ",";
                cout << "x" << j;
            }
            cout << "]) / (x" << i << " - x0)" << endl;

            double numerator = table[0][i];
            cout << "      = (" << table[1][i - 1] << " - " << table[0][i - 1] << ") / (" << x[i] << " - " << x[0] << ") = "
                << numerator << endl;
        }
    }

    cout << "\n Сборка полинома Ньютона шаг за шагом:\n";
    cout << "P(x) = ";

    for (int i = 0; i < n; i++)
    {
        if (i > 0) cout << "\n     + ";
        cout << table[0][i];

        for (int j = 0; j < i; j++)
        {
            cout << " * (x - " << x[j] << ")";
        }

        if (i == 0)
            cout << "     <--  Первый член: просто y при x = x0";
        else
            cout << "     <--  Член " << i << ": включает множители (x - xj)";
    }

    cout << "\n\n Полином в числовом виде:\nP(x) = ";

    for (int i = 0; i < n; i++)
    {
        if (i > 0) cout << " + ";
        cout << table[0][i];
        for (int j = 0; j < i; j++)
        {
            cout << "*(x - " << x[j] << ")";
        }
    }

    cout << "\n";
}

// Пошаговое вычисление значения полинома
void rashet_polinoma_v_tochke(const vector<double>& x, const vector<vector<double>>& table, double value) {
    int n = x.size();
    double result = table[0][0];
    double term = 1.0;

    cout << "\n Пошаговое вычисление значения полинома в точке x = " << value << ":\n";
    cout << "Шаг 0: Берём первый коэффициент: " << table[0][0] << "\n";

    for (int i = 1; i < n; i++)
    {
        term = 1.0;
        cout << "\nШаг " << i << ": Вычисляем множитель: ";

        for (int j = 0; j < i; j++)
        {
            cout << "(x - x[" << j << "]) = (" << value << " - " << x[j] << ")";
            if (j < i - 1) cout << " * ";
        }

        for (int j = 0; j < i; j++)
        {
            term *= (value - x[j]);
        }

        double add = term * table[0][i];
        cout << "\n  Произведение: " << term;
        cout << "\n  Умножаем на коэффициент: " << table[0][i];
        cout << "\n  Прибавляем к результату: " << term << " * " << table[0][i] << " = " << add << "\n";

        result += add;
    }

    cout << "\n Окончательный ответ: P(" << value << ") = " << result << "\n";
}

int main()
{
    setlocale(LC_ALL, "Russian");

    int n;
    string input_n;
    cout << "Введите количество узлов интерполяции: ";
    getline(cin, input_n);
    stringstream ss(input_n);

    if (!(ss >> n) || !(ss.eof()) || n <= 0)
    {
        cerr << "Ошибка: введите положительное целое число без лишних символов.\n";
        return 1;
    }
    if (n > MAX_YZEL)
    {
        cerr << "Ошибка: превышено максимально допустимое количество узлов (" << MAX_YZEL << ").\n";
        return 1;
    }
    if (n > 20)
    {
        cout << "Предупреждение: при n > 20 возможно снижение точности и увеличение времени вычислений.\n";
    }

    vector<double> x(n), y(n);
    cout << "\nВведите пары значений x и y построчно:\n";
    for (int i = 0; i < n; i++)
    {
        cout << "x[" << i << "] и y[" << i << "]: ";
        cin >> x[i] >> y[i];

        if (cin.fail())
        {
            cerr << "Ошибка: некорректный ввод. Ожидаются числовые значения.\n";
            return 1;
        }

        for (int j = 0; j < i; j++)
        {
            if (abs(x[i] - x[j]) < MIN_DIFF)
            {
                cerr << "Ошибка: значения x[" << i << "] и x[" << j << "] слишком близки или совпадают.\n";
                return 1;
            }
        }
    }

    cout << "\n Введённые данные:\n";
    cout << setw(10) << "x" << setw(10) << "y" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << setw(10) << x[i] << setw(10) << y[i] << endl;
    }

    vector<vector<double>> table = tabl_razdelen_razn(x, y);
    printPolinom(x, table);

    double value;
    cout << "\nВведите точку, в которой нужно вычислить значение полинома: ";
    cin >> value;

    if (cin.fail())
    {
        cerr << "Ошибка: некорректный ввод значения x.\n";
        return 1;
    }

    rashet_polinoma_v_tochke(x, table, value);

    return 0;
}
