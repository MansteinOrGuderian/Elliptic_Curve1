#include "Header.h"

unsigned int mod(long value, unsigned mod) {
	while (value < 0)
		value += mod;
	value %= mod;
	return static_cast<unsigned int>(value);
}

unsigned int inverse(unsigned &number_to_inverse, unsigned modulo_number) { //  number^-1 mod mod
	if (number_to_inverse == 0 || modulo_number == 0)// gcd(0, 0) = 0, gcd (0, n) = n, gcd (n, 0) = n
		return 0; // In this cases inverse don't exist
	if (number_to_inverse > modulo_number)
		number_to_inverse = number_to_inverse % modulo_number;
	long long int r_i_minus_one = static_cast<long long>(number_to_inverse), r_i = static_cast<long long>(modulo_number), q(0), r_i_plus_one(1);
	std::vector<long long int> r_i_minus_one_vector, r_i_vector, q_vector, minus_q_vector, r_i_plus_one_vector;
	while (r_i_plus_one != 0) { // common Euclidean Algorithm
		q = r_i_minus_one / r_i;
		r_i_plus_one = r_i_minus_one % r_i;
		q_vector.push_back(q), minus_q_vector.push_back(-q), r_i_plus_one_vector.push_back(r_i_plus_one);
		r_i_minus_one_vector.push_back(r_i_minus_one), r_i_vector.push_back(r_i);
		r_i_minus_one = r_i;
		r_i = r_i_plus_one;
	}
	//for (long long int i = 0; i < q_vector.size(); i++)
	//	std::cout << r_i_minus_one_vector[i] << " = " << r_i_vector[i] << " * " << q_vector[i] << " + " << r_i_plus_one_vector[i] << '\n';
	if (r_i_vector.back() != 1) { // gcd(a, b) != 1, a^{-1} don't exist
		//std::cout << "Doesn't exict a^{-1} mod b\n" << number_to_inverse << "^{-1} mod " << modulo_number << '\n';
		return 0; // error code
	}
	else {
		//std::cout << "Exict a^{-1} mod b\n" << number_to_inverse << "^{-1} mod " << modulo_number << '\n';
		minus_q_vector.erase(minus_q_vector.begin()); // first element in vector always 0, so we delete it
		long long int temp_i_minus_one = 0, temp_i = 1, temp_i_plus_one;
		std::vector<long long int> temp_vector;
		temp_vector.push_back(0), temp_vector.push_back(1);
		for (long long int i = 0; i < minus_q_vector.size(); i++) {
			temp_i_plus_one = minus_q_vector[i] * temp_i + temp_i_minus_one;
			temp_i_minus_one = temp_i;
			temp_i = temp_i_plus_one;
			temp_vector.push_back(temp_i_plus_one);
		}
		//for (long long int i = 0; i < minus_q_vector.size(); i++) {
		//	if (i == 0)
		//		std::cout << " |  |";
		//	std::cout << minus_q_vector[i] << ((i != minus_q_vector.size() - 1) ? '|' : '\n');
		//}
		//for (long long int i = 0; i < temp_vector.size(); i++)
		//	std::cout << temp_vector[i] << ((i != temp_vector.size() - 1) ? '|' : '\n');
		unsigned inversed = (temp_vector[(temp_vector.size() - 2)] < 0 ? static_cast<unsigned>(temp_vector[(temp_vector.size() - 2)] + static_cast<long long>(modulo_number)) : static_cast<unsigned>(temp_vector[(temp_vector.size() - 2)]));
		return inversed; // need penultimate element, because last element is modulo
	}
}

template <typename T>
T power(T value, T pow, T mod) {
	T temp_value = static_cast<T>(1);
	for (size_t i = 0; i < pow; i++) {
		temp_value *= value;
		temp_value %= mod;
	}
	return temp_value;
}

std::ostream& operator<<(std::ostream& os, const Point& cur_point) {
	os << '(' << cur_point.x << ", " << cur_point.y << ')';
	return os;
}

Point Elliptic_Curve::sum_neq_points(Point P, Point Q) { // P (x_1, y_1), Q (x_2, y_2)
	Point result; // default (0, 0)
	if (P.x != Q.x) {
		unsigned int temp_numerator = mod(static_cast<long>(P.y - Q.y), PRIME_P);
		unsigned int temp_denumerator = mod(static_cast<long>(P.x - Q.x), PRIME_P);
		temp_denumerator = inverse(temp_denumerator, PRIME_P);
		unsigned temp = mod(temp_numerator * temp_denumerator, PRIME_P);

		// x = [(y_1 - y_2)/(x_1 - x_2)]^2 - x_1 - x_2
		result.x = mod(static_cast<long>(mod(static_cast<long>(mod(power(temp, static_cast<unsigned>(2), PRIME_P), PRIME_P)) - P.x, PRIME_P)) - Q.x, PRIME_P); // Cringe

		// y = [(y_1 - y_2)/(x_1 - x_2)] * (x_1 - x_3) - y_1
		result.y = mod(static_cast<long>(mod(temp * mod(static_cast<long>(P.x - result.x), PRIME_P), PRIME_P) - P.y), PRIME_P); // Cringe
		return result;
	}
	else // P.x = Q.x -> it's inverse points. P = -Q
		return result; //double_Point(P); // infinity point
}

Point Elliptic_Curve::double_Point(Point P) { // P (x_1, y_1)
	Point result;
	if (P.y == 0)
		return result; // infinity point 
	unsigned int temp_numerator = mod(mod(3 * mod(power(P.x, static_cast<unsigned>(2), PRIME_P), PRIME_P), PRIME_P) + Elliptic_Curve::a, PRIME_P);
	unsigned int temp_denumerator = mod(2 * P.y, PRIME_P);
	temp_denumerator = inverse(temp_denumerator, PRIME_P);
	unsigned temp = mod(temp_numerator * temp_denumerator, PRIME_P);

	// x = [(3 x_1^2 + a)/(2 * y_1)]^2 - 2 * x_1
	result.x = mod(mod(power(temp, static_cast<unsigned>(2), PRIME_P), PRIME_P) - mod(2 * P.x, PRIME_P), PRIME_P); // Cringe

	// y = [(3 x_1^2 + a)/(2 * y_1)] * (x_1 - x_3) - y_1
	result.y = mod(static_cast<long>(mod(temp * mod(static_cast<long>(P.x - result.x), PRIME_P), PRIME_P) - P.y), PRIME_P); // Cringe

	return result;
} // P = Q (2P)

Point Elliptic_Curve::minus_Point(Point P) {
	Point result(P.x, mod(static_cast<long>((-1) * P.y), PRIME_P));
	return result; // x = x;  y = (-y) mod mod
} // Q = -P


void Elliptic_Curve::get_all_curve_points() { // y^2 = x^3 + ax + b
	std::set<unsigned> possible_quad_residue;
	unsigned *temp_y2_array = new unsigned [PRIME_P] ();
	for (unsigned cur_pos = 0; cur_pos < PRIME_P; cur_pos++) {
		if (cur_pos > 0) {
			possible_quad_residue.insert(mod(power(cur_pos, static_cast<unsigned>(2), PRIME_P), PRIME_P)); // Could be optimised, cause after (PRIME_P - 1) in reverse order
			//std::cout << mod(power(cur_pos, static_cast<unsigned>(2)), PRIME_P) << ' ';
		}
		temp_y2_array[cur_pos] = mod(mod(mod(power(cur_pos, static_cast<unsigned>(3), PRIME_P), PRIME_P) + mod(a * cur_pos, PRIME_P), PRIME_P) + b, PRIME_P);
	}

	std::cout << "All possible quadratic residues:\n";
	for (auto element : possible_quad_residue) {
		std::cout << element << " ";
	}
	std::cout << '\n';

	bool* is_quad_residue = new bool[PRIME_P]();

	for (unsigned cur_pos = 0; cur_pos < PRIME_P; cur_pos++)
		if(possible_quad_residue.find(temp_y2_array[cur_pos]) != possible_quad_residue.end())
			is_quad_residue[cur_pos] = true;

	points_collection.push_back(Point()); // (0, 0) - infinity point

	for (unsigned cur_pos = 0; cur_pos < PRIME_P; cur_pos++) { // works ONLY for p = 4k+3 (solving y^2 = a mod p)
		if (is_quad_residue[cur_pos]) {
			unsigned k = (PRIME_P - 3) / 4; 
			unsigned y1 = mod(power(temp_y2_array[cur_pos], k + static_cast<unsigned>(1), PRIME_P), PRIME_P); // y1 = a^(k+1) mod p
			unsigned y2 = mod(static_cast<long>(-1 * y1), PRIME_P); // y2 = -y1 mod p

			Point first(cur_pos, y1), second(cur_pos, y2);
			points_collection.push_back(first);
			points_collection.push_back(second);
		}
		if (temp_y2_array[cur_pos] == 0)
			points_collection.push_back(Point(cur_pos, static_cast<unsigned>(0)));
	}

	// --- START TABLE PART ---
	const unsigned COLUMNS_PER_CHUNK = 15;
	const unsigned LABEL_WIDTH = 7; // Width for "| x |" or "| RHS |" label column

	std::cout << "\n--- Elliptic Curve Points Table ---\n";
	std::cout << "Curve: y^2 = x^3 + " << a << "x + " << b << " (mod " << PRIME_P << ")\n\n";

	// Loop through the data in chunks of COLUMNS_PER_CHUNK
	for (unsigned start_x = 0; start_x < PRIME_P; start_x += COLUMNS_PER_CHUNK) {

		// Determine the end of the current chunk
		unsigned end_x = start_x + COLUMNS_PER_CHUNK;
		if (end_x > PRIME_P) {
			end_x = PRIME_P;
		}

		// Calculate the total width needed for the separator line
		// Label width (7) + Inner Pipe (|) (1) + (8 chars per data column * N columns)
		// 8 chars = setw(6) + space (1) + pipe (|) (1)
		unsigned separator_length = LABEL_WIDTH + 1 + (8 * (end_x - start_x));
		// Function to print a separator line of correct width
		auto print_separator = [&]() {
			// Print the starting pipe '|'
			std::cout << "|";
			// Print hyphens for the exact width of the content (separator_length - 1 characters)
			for (unsigned i = 0; i < separator_length - 1; ++i) {
				std::cout << "-";
			}
			// Print the closing pipe '|'
			std::cout << "|\n";
		};

		// --- Column Headers (x values) ---
		print_separator();
		std::cout << "|" << std::setw(LABEL_WIDTH) << "x" << "|";
		for (unsigned x = start_x; x < end_x; x++) {
			std::cout << std::setw(6) << x << " |";
		}
		std::cout << "\n";
		print_separator();

		// --- Row 1: x^3 + ax + b (RHS value) ---
		std::cout << "|" << std::setw(LABEL_WIDTH) << "RHS" << "|";
		for (unsigned x = start_x; x < end_x; x++) {
			std::cout << std::setw(6) << temp_y2_array[x] << " |";
		}
		std::cout << "\n";
		print_separator();

		// --- Row 2: y1 (first square root) ---
		std::cout << "|" << std::setw(LABEL_WIDTH) << "y1" << "|";
		for (unsigned x = start_x; x < end_x; x++) {
			if (is_quad_residue[x] || temp_y2_array[x] == 0) {
				unsigned rhs_value = temp_y2_array[x];
				unsigned y_val;

				if (rhs_value == 0) {
					y_val = 0;
				}
				else {
					unsigned k = (PRIME_P - 3) / 4;
					y_val = mod(power(rhs_value, k + 1, PRIME_P), PRIME_P);
				}
				std::cout << std::setw(6) << y_val << " |";
			}
			else {
				std::cout << std::setw(6) << "-" << " |";
			}
		}
		std::cout << "\n";
		print_separator();

		// --- Row 3: y2 (second square root) ---
		std::cout << "|" << std::setw(LABEL_WIDTH) << "y2" << "|";
		for (unsigned x = start_x; x < end_x; x++) {
			if (is_quad_residue[x]) {
				unsigned rhs_value = temp_y2_array[x];
				unsigned y_val;

				if (rhs_value == 0) {
					y_val = 0;
				}
				else {
					unsigned k = (PRIME_P - 3) / 4;
					unsigned y1 = mod(power(rhs_value, k + 1, PRIME_P), PRIME_P);
					y_val = mod(static_cast<long>(-1 * y1), PRIME_P);
				}
				std::cout << std::setw(6) << y_val << " |";
			}
			else {
				std::cout << std::setw(6) << "-" << " |";
			}
		}
		std::cout << "\n";
		print_separator();

		std::cout << "\n"; // Add a gap between chunks for better readability
	}
}

std::vector<bool> decimal_to_binary(unsigned int decimal_number) {
	std::vector<bool> binary_vector;
	if (decimal_number == 0)
		binary_vector.push_back(false);
	while (decimal_number > 0) {
		// Check the least significant bit (LSB) using the bitwise AND operator (&)
		// If the LSB is 1, decimal_number & 1 == 1 (true), otherwise 0 (false).
		binary_vector.push_back(decimal_number & 1);
		decimal_number >>= 1;
	}
	return binary_vector; // start from 2^0 (least significant bit (LSB))
}

Point Elliptic_Curve::straight_line(Point P, Point Q) {
	Point pair; 
	//std::cout << "\n(x - x_Q) / (x_P - x_Q) = (y - y_Q) / (y_P - y_Q)\n";
	unsigned denominator_one = mod(static_cast<long>(P.x - Q.x), PRIME_P);
	unsigned inverse_denominator_one = inverse(denominator_one, PRIME_P);
	unsigned denominator_two = mod(static_cast<long>(P.y - Q.y), PRIME_P);
	pair.x = mod(static_cast<long>(inverse_denominator_one*denominator_two), PRIME_P); 
	pair.y = mod(mod(static_cast<long>(pair.x) * (-1) * static_cast<long>(Q.x), PRIME_P) + Q.y, PRIME_P); 
	// std::cout << "Straight:\ny = " << pair.x << "*x + " << pair.y << '\n';
	return pair; // k == pair.x; b == pair.y;    y = kx + b
}

Point Elliptic_Curve::tangent_line(Point P) {
	Point pair;
	//std::cout << "\n(3 * x_P ^2 + curve.a) / (2 y_P)  * (x - x_P) + y_P\n";
	unsigned denominator = mod(static_cast<long>(2 * P.y), PRIME_P);
	unsigned inverse_denominator = inverse(denominator, PRIME_P);
	pair.x = mod(static_cast<long>(inverse_denominator * mod(3 * power(P.x, static_cast<unsigned>(2), PRIME_P) + a, PRIME_P)), PRIME_P);
	pair.y = mod(mod(static_cast<long>(pair.x) * (-1) * static_cast<long>(P.x), PRIME_P) + P.y, PRIME_P);
	// std::cout << "Tangent:\ny = " << pair.x << "*x + " << pair.y << '\n';
	return pair; // k == pair.x; b == pair.y;    y = kx + b
}

Point Elliptic_Curve::vertical_line(Point P) {
	Point pair (P.x, static_cast<unsigned>(0));
	// std::cout << "Vertical:\nx = " << P.x << '\n';
	return pair; // x = x_P, y = 0
}

unsigned int Elliptic_Curve::Millers_algorithm(Point P) {
	Point T = P;
	unsigned f = 1;
	std::vector<bool> binary_prime = decimal_to_binary(static_cast<unsigned>(5));
	for (int j = binary_prime.size() - 1; j >= 0; j--) {
		Point I_T_T = tangent_line(T);
		Point v_2T = vertical_line(double_Point(T));
		unsigned fraction = mod(static_cast<long>(mod(static_cast<long>(I_T_T.x * T.x), PRIME_P) + I_T_T.y), PRIME_P);
		fraction = mod(static_cast<long>(fraction * inverse(v_2T.x, PRIME_P)), PRIME_P);
		f = mod(static_cast<long>(power(f, static_cast<unsigned>(2), PRIME_P) * fraction), PRIME_P);
		T = double_Point(T);
		if (binary_prime[j] == true) {
			Point I_T_P = straight_line(T, P);
			Point v_T_plus_P = vertical_line(sum_neq_points(T, P));
			unsigned fraction_temp = mod(static_cast<long>(mod(static_cast<long>(I_T_P.x * T.x), PRIME_P) + I_T_P.y), PRIME_P);
			fraction_temp = mod(static_cast<long>(fraction_temp * inverse(v_T_plus_P.x, PRIME_P)), PRIME_P);
			f = mod(static_cast<long>(f * fraction_temp), PRIME_P);
			T = sum_neq_points(T, P);
		}
	}
	//std::vector<bool> binary = decimal_to_binary(a);
	//for (bool bit : binary) {
	//	std::cout << (bit ? '1' : '0');
	//}
	return f;
}


int main() {
	Elliptic_Curve x2_30x_34;
	//x2_30x_34.get_all_curve_points();
	//std::cout << "In total " << x2_30x_34.points_collection.size() << " points in Elliptic Curve. Here are they:\n";
	//for (int i = 0; i < x2_30x_34.points_collection.size(); i++)
	//	std::cout << x2_30x_34.points_collection[i] << '\n';// << (i % 10 == 0 ? '\n' : ' ');
	//
	//std::cout << '\n';
	std::vector<Point> Ps_vector;
	Point P(36, 571);
	std::cout << "P " << P << '\n';
	Ps_vector.push_back(P);
	Point P2 = x2_30x_34.double_Point(P);
	std::cout << "2P " << P2 << '\n';
	Ps_vector.push_back(P2);
	Point P3 = x2_30x_34.sum_neq_points(P, P2);
	std::cout << "3P " << P3 << '\n';
	Ps_vector.push_back(P3);
	Point P4 = x2_30x_34.double_Point(P2);
	std::cout << "4P " << P4 << '\n';
	Ps_vector.push_back(P4);
	Point P5 = x2_30x_34.sum_neq_points(P, P4);
	std::cout << "5P " << P5 << '\n';
	Ps_vector.push_back(P5);
	std::cout << '\n';

	/*Point P6 = x2_30x_34.double_Point(P3);
	std::cout << "6P " << P6 << '\n';
	Point P7 = x2_30x_34.sum_neq_points(P3, P4);
	std::cout << "7P " << P7 << '\n';
	Point P8 = x2_30x_34.double_Point(P4);
	std::cout << "8P " << P8 << '\n';
	Point P9 = x2_30x_34.sum_neq_points(P3, P6);
	std::cout << "9P " << P9 << '\n';
	Point P10 = x2_30x_34.sum_neq_points(P2, P8);
	std::cout << "10P " << P10 << '\n';
	std::cout << '\n';*/

	std::vector<Point> Qs_vector;
	Point Q(420, 48); // Q (121, 244)
	std::cout << "Q " << Q << '\n';
	Qs_vector.push_back(Q);
	Point Q2 = x2_30x_34.double_Point(Q);
	std::cout << "2Q " << Q2 << '\n';
	Qs_vector.push_back(Q2);
	Point Q3 = x2_30x_34.sum_neq_points(Q, Q2);
	std::cout << "3Q " << Q3 << '\n';
	Qs_vector.push_back(Q3);
	Point Q4 = x2_30x_34.double_Point(Q2);
	std::cout << "4Q " << Q4 << '\n';
	Qs_vector.push_back(Q4);
	Point Q5 = x2_30x_34.sum_neq_points(Q, Q4);
	std::cout << "5Q " << Q5 << '\n';
	//Qs_vector.push_back(Q5);
	std::cout << '\n';

	std::vector<Point> P_Q_vector;
	int counter = 0;
	for (int p_i = 0; p_i < 4; p_i++) {
		P_Q_vector.push_back(Ps_vector[p_i]);
		for (int q_j = 0; q_j < 4; q_j++) {
			if (p_i == 0) {
				P_Q_vector.push_back(Ps_vector[q_j]);
			}
			std::cout << ((p_i + 1) == 1 ? ' ' : char('0' + p_i + 1)) << "P + " << ((q_j + 1) == 1 ? ' ' : char('0' + q_j + 1)) << "Q = " << Ps_vector[p_i] << " + " << Qs_vector[q_j] << " = ";
			Point result = x2_30x_34.sum_neq_points(Ps_vector[p_i], Qs_vector[q_j]);
			std::cout << result << '\n';
			P_Q_vector.push_back(result);
			counter++;
		}
	}
	std::cout << "Amount of points in subgroup E[5] exponent 5 is: " << Ps_vector.size() + Qs_vector.size() + counter << '\n';

	Point R(0, 36);		
	std::vector<Point> R_current;
	R_current.push_back(R);
	std::cout << "R = " << R_current[0] << "\n";
	for (int k = 2; k <= 130; ++k) {
		Point R_k_minus_1 = R_current.back(); // P_k_minus_1 is (k-1)R
		Point R_k;

		if (R_k_minus_1.x == R.x && R_k_minus_1.y == R.y) {
			// Avoids calling sum_neq_points when adding two identical points.
			// This happens when k-1 = 1, so P_k_minus_1 = R, and we calculate 2R.
			R_k = x2_30x_34.double_Point(R);
		}
		else {
			// General case: kR = (k-1)R + R
			R_k = x2_30x_34.sum_neq_points(R_k_minus_1, R);
		}

		// Store the result for the next iteration
		R_current.push_back(R_k);

		// --- Output kR ---
		std::cout << std::to_string(k) << "R = "
			<< std::to_string(k - 1) << "R + "
			<< "R" // The second point is always R (1R)
			<< " = " << R_k << "\n";
	}

	std::cout << "Is subgroup R and P_Q intersect check:\n";
	for (auto P_Q : P_Q_vector)
		if (std::find(R_current.begin(), R_current.end(), P_Q) != R_current.end())
			std::cout << P_Q << " Same found\n";
	std::cout << "Groups do not intersect\n";

	unsigned Q_S = x2_30x_34.Millers_algorithm(x2_30x_34.sum_neq_points(Q, R));
	std::cout << Q_S << ' ';
	unsigned S = x2_30x_34.Millers_algorithm(R);
	std::cout << S << ' ';
	unsigned Weil_pairing = mod(static_cast<long>(Q_S * inverse(S, static_cast<unsigned>(631))), static_cast<unsigned>(631));
	std::cout << Weil_pairing << ' ' << power(Weil_pairing, static_cast<unsigned>(5), static_cast<unsigned>(631));


	//Point Q6 = x2_30x_34.double_Point(Q3);
	//std::cout << "6Q " << Q6 << '\n';
	//Point Q7 = x2_30x_34.sum_neq_points(Q3, Q4);
	//std::cout << "7Q " << Q7 << '\n';
	//Point Q8 = x2_30x_34.double_Point(Q4);
	//std::cout << "8Q " << Q8 << '\n';
	//Point Q9 = x2_30x_34.sum_neq_points(Q3, Q6);
	//std::cout << "9Q " << Q9 << '\n';
	//Point Q10 = x2_30x_34.sum_neq_points(Q3, Q7);
	//std::cout << "10Q " << Q10 << '\n';
	//std::cout << '\n';

	//Point sum = x2_30x_34.sum_neq_points(P, Q);
	//std::cout << sum;

}
