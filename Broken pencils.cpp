
#include <iostream>
#include <vector>
#include <fstream>

#define pi 3.14159265359

const int n = 100;				// Кількість ітерацій для інтегрування
const double gamma1 = 1.0E-9;	// Поверхнева енергія каоліну (НЕВІДОМА!)
const double gamma2 = 1.0E-9;	// Поверхнева енергія графіту (НЕВІДОМА!)
const double E1 = 2.5E9;		// Модуль Юнга каоліну = 2.5 ГН/м2 (при Т = 273 К)
const double E2 = 5.88E9;		// Модуль Юнга графіта = 5880 МН/м2 (при Т = 293 К)
const double alfa = 2;			// Коефіцієнт в формулі Гріффітса
const double ALFA = 0.15;		// Коефіцієнт конвертації потенціальної енергії падіння в поверхневу енергію

const double D = 0.018;			// Довжина грифеля в поперечному перерізі
const double D2 = 0.026;
const double L = 0.176;
const double ro0 = 1.3;
const double m = 3.5E-3;		// Маса олівця
const double Cx = 1.0;

const int width = 10;			// Ширина модельованої області (ширина грифеля, виражена в кількості вузлів)
const int height = 36;			// Висота модельованої області
const int center_X = 5;			// Х-координата центру руйнування
const int center_Y = height / 2;		// Y-координата центру руйнування

double rat = 0.5;
double length_kernel = 1.*D / (width * (1 + rat));	// Довжина графітового зерна
double length_bonding = 1.*rat * length_kernel;		// Відстань між зернами графіту
double step_bonding = 1.*length_bonding / n;
double step_kernel = 1.*length_kernel / n;


/**
Кількість шляхів деякої чверті координатної площини
good_ways1 означає хороші шляхи 1 чверті.
Хороші шляхи з'єднують початкову точку (центр руйнування) з лівим або правим краєм області.
Погані шляхи - всі інші.
*/
int good_ways1 = 0;
int bad_ways1 = 0;
int good_ways2 = 0;
int bad_ways2 = 0;
int good_ways3 = 0;
int bad_ways3 = 0;
int good_ways4 = 0;
int bad_ways4 = 0;

std::vector<int> way_length1;
std::vector<int> way_length2;
std::vector<int> way_length3;
std::vector<int> way_length4;

/**
Ця функція запускає алгоритм проходу по графу, починаючи від точки з координатами (center_X; center_Y)
тобто центру руйнування. Проходить алгоритм у 4 напрямках відповідно до чверті на координатній площині.
В результаті отримується число шляхів в кожній чверті і довжина кожного шляху.
*/
int searchWays(int x, int y, int mode) {
	if (x == width && center_X != width) {
		int l = (width - center_X) + abs(y - center_Y);
		if (y >= center_Y) {
			good_ways1++;
			way_length1.push_back(l);
		}
		else {
			good_ways4++;
			way_length4.push_back(l);
		}
		return 0;
	}
	if (x == 0 && center_X != 0) {
		int l = center_X + abs(y - center_Y);
		if (y >= center_Y) {
			good_ways2++;
			way_length2.push_back(l);
		}
		else {
			good_ways3++;
			way_length3.push_back(l);
		}
		return 0;
	}
	if (y == height && center_Y != height) {
		if (x >= center_X) bad_ways1++;
		else bad_ways2++;
		return 0;
	}
	if (y == 0 && center_Y != 0) {
		if (x >= center_X) bad_ways4++;
		else bad_ways3++;
		return 0;
	}

	switch (mode) {
	case(1):
		searchWays(x + 1, y, 1);
		searchWays(x, y + 1, 1);
		break;
	case(2):
		searchWays(x - 1, y, 2);
		searchWays(x, y + 1, 2);
		break;
	case(3):
		searchWays(x - 1, y, 3);
		searchWays(x, y - 1, 3);
		break;
	case(4):
		searchWays(x + 1, y, 4);
		searchWays(x, y - 1, 4);
		break;
	default:
		std::cout << "Invalid mode";
	}
}

/**
Ця функція розраховує енергію, яку потрібно прикласти, щоб розірвати nodes вузлів.
Для розрахунку функція використовує формулу Гріффітса, інтегрує критичну напругу
від значення довжини lendth до length + bond_length а потім від
(length + bond_length) до (length + bond_length + kernel_length),
множить на поперечний переріз вузла (kernel_length * pi/4).
Ця операція застосовна до кожного окремого вузла, а значення length є біжучим.
*/
double calculateStrain(std::vector<int> way) {

	double l = length_kernel;
	double function = 0.0;

	double totalStrain = 0.0;

	double x = center_X * (length_bonding + length_kernel);
	for (int k = 0; k < way.size(); k++) {
		x += way[k] * (length_bonding + length_kernel);

		double initial_f = sqrt(gamma1 * E1 / (alfa * (l)));
		for (double i = step_bonding; i < length_bonding; i += step_bonding) {
			function = sqrt(gamma1 * E1 / (alfa * (l + i)));
			totalStrain += step_bonding * sqrt(x * (D + 1E-6) - x * x)  * (function + initial_f);
			initial_f = function;
		}

		l += length_bonding;
		initial_f = sqrt(gamma2 * E2 / (alfa * (l)));
		for (double i = step_kernel; i < length_kernel; i += step_kernel) {
			function = sqrt(gamma2 * E2 / (alfa * (l + i)));
			totalStrain += step_kernel * sqrt(x*(D + 1E-6) - x * x) * (function + initial_f);
			initial_f = function;
		}
		l += length_kernel;
	}

	return totalStrain;
}

int sum(std::vector<int> way) {
	int acc = 0;
	for (int i : way) acc += i;
	return acc;
}

int zeros(std::vector<int> way) {
	int count = 0;
	for (int i : way) {
		if (i == 0) count++;
	}
	return count;
}

/**
Ця функція розраховує, які із знайдених шляхів можуть бути реалізовані енергетично.
Для цього визначається яка максимальна кількість вузлів може бути пройдена
і підраховується кількість шляхів, для яких n <= nodes
Результатом цієї функції є кількість енергетично можливих шляхів.
*/
int numberOfSuitableWays(double Energy) {
	int number_of_suitable_ways = 0;
	double applied_Energy = ALFA * Energy;

	int y = 0;
	while (y < height - center_Y) {
		std::vector<int> way;
		int local_y = y;
		way.push_back(1);
		for (int j = 1; j < (width - center_X) + y; j++) {
			if (j > (width - y) / 2 - center_X && local_y > 0) {
				way.push_back(0);
				local_y--;
			}
			else way.push_back(1);
		}
		double nesessary_Energy = calculateStrain(way);
		if (nesessary_Energy > applied_Energy) {
			break;
		}
		y++;
	}
	for (int a : way_length1) {
		if (a > (width + 2 * y)) continue;
		for (int b : way_length3) {
			if (a + b <= (width + 2 * y)) number_of_suitable_ways++;
		}
	}
	for (int a : way_length2) {
		if (a > (width + 2 * y)) continue;
		for (int b : way_length4) {
			if (a + b <= (width + 2 * y)) number_of_suitable_ways++;
		}
	}
	return number_of_suitable_ways;
}

int main() {
	searchWays(center_X, center_Y, 1);
	searchWays(center_X, center_Y, 2);
	searchWays(center_X, center_Y, 3);
	searchWays(center_X, center_Y, 4);

	double good_ways = good_ways1 * good_ways3 + good_ways2 * good_ways4;
	double bad_ways = bad_ways1 * bad_ways3 + bad_ways2 * bad_ways4;
	double all_ways = good_ways + bad_ways;
	double ratio = 1. * good_ways / all_ways;
	std::cout << "all_ways = " << all_ways << std::endl;
	std::cout << "good_ways = " << good_ways << std::endl;
	std::cout << "ratio = " << ratio << std::endl;

	std::ofstream fout;
	fout.open("Table.txt");

	std::cout << "Height (m) \t Probability (%)\n\n";
	for (double h = 0; h < 3; h += 0.05) {
		double beta = Cx * (D2 * L) * ro0 / 2;
		double V = -beta * h / m + sqrt(pow(beta*h / m, 2) + 2 * 9.8*h);
		double Energy = 0.5*m*V*V;
		double number_of_suitable_ways = numberOfSuitableWays(Energy);
		std::cout << h << " \t " << V << " \t " << 100.*number_of_suitable_ways / good_ways << std::endl;
		fout << h << " \t " << 100.*number_of_suitable_ways / good_ways * ratio << std::endl;
	}

	system("pause");
	return 0;
}
