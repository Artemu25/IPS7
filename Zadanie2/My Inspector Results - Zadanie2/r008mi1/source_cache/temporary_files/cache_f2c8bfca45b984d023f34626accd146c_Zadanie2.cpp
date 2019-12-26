#include "pch.h"
#include <vector>
#include <thread>
#include <stdio.h>
#include <exception>
#include <locale.h>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>


#include <iostream>

/// ������������, ������������ ��� ����� ����������� ����������
/// ������� �������� �������: �� ������� ��� �� ��������
enum class eprocess_type
{
	by_rows = 0,
	by_cols
};

void InitMatrix(double** matrix, const size_t numb_rows, const size_t numb_cols)
{
	for (size_t i = 0; i < numb_rows; ++i)
	{
		for (size_t j = 0; j < numb_cols; ++j)
		{
			matrix[i][j] = rand() % 5 + 1;
		}
	}
}

/// ������� PrintMatrix() �������� �������� ������� <i>matrix</i> �� �������;
/// numb_rows - ���������� ����� � �������� ������� <i>matrix</i>
/// numb_cols - ���������� �������� � �������� ������� <i>matrix</i>
void PrintMatrix(double** matrix, const size_t numb_rows, const size_t numb_cols)
{
	printf("Generated matrix:\n");
	for (size_t i = 0; i < numb_rows; ++i)
	{
		for (size_t j = 0; j < numb_cols; ++j)
		{
			printf("%lf ", matrix[i][j]);
		}
		printf("\n");
	}
}

/// ������� FindAverageValues() ������� ������� �������� � ������� <i>matrix</i>
/// �� �������, ���� �� �������� � ����������� �� �������� ��������� <i>proc_type</i>;
/// proc_type - �������, � ����������� �� �������� ������� �������� ����������� 
/// ���� �� �������, ���� �� ������� �������� ������� <i>matrix</i>
/// matrix - �������� �������
/// numb_rows - ���������� ����� � �������� ������� <i>matrix</i>
/// numb_cols - ���������� �������� � �������� ������� <i>matrix</i>
/// average_vals - ������, ���� ����������� ����������� ������� ��������
void FindAverageValues(eprocess_type proc_type, double** matrix, const size_t numb_rows, const size_t numb_cols, double* average_vals)
{
	switch (proc_type)
	{
	case eprocess_type::by_rows:
	{
		for (size_t i = 0; i < numb_rows; ++i)
		{
			double sum(0.0);
			for (size_t j = 0; j < numb_cols; ++j)
			{
				sum += matrix[i][j];
			}
			average_vals[i] = sum / numb_cols;
		}
		break;
	}
	case eprocess_type::by_cols:
	{
		for (size_t j = 0; j < numb_cols; ++j)
		{
			double sum(0.0);
			for (size_t i = 0; i < numb_rows; ++i)
			{
				sum += matrix[i][j];
			}
			average_vals[j] = sum / numb_rows;
		}
		break;
	}
	default:
	{
		throw("Incorrect value for parameter 'proc_type' in function FindAverageValues() call!");
	}
	}
}


void FindAverageValuesCilk(eprocess_type proc_type, double** matrix, const size_t numb_rows, const size_t numb_cols, double* average_vals)
{
	switch (proc_type)
	{
	case eprocess_type::by_rows:
	{
		cilk_for(size_t i = 0; i < numb_rows; ++i)
		{
			double sum(0.0);
			for (size_t j = 0; j < numb_cols; ++j)
			{
				sum += matrix[i][j];
			}
			average_vals[i] = sum / numb_cols;
		}
		break;
	}
	case eprocess_type::by_cols:
	{
		cilk_for(size_t j = 0; j < numb_cols; ++j)
		{
			double sum(0.0);
			for (size_t i = 0; i < numb_rows; ++i)
			{
				sum += matrix[i][j];
			}
			average_vals[j] = sum / numb_rows;
		}
		break;
	}
	default:
	{
		throw("Incorrect value for parameter 'proc_type' in function FindAverageValues() call!");
	}
	}
}

/// ������� PrintAverageVals() �������� �������� ������� <i>average_vals</i> �� �������;
/// proc_type - �������, ���������� �� ��, ��� ���� ��������� 
/// ������� �������� �������� ������� �� ������� ��� �� ��������
/// average_vals - ������, �������� ������� �������� �������� �������,
/// ����������� �� ������� ��� �� ��������
/// dimension - ���������� ��������� � �������� ������� <i>average_vals</i>
void PrintAverageVals(eprocess_type proc_type, double* average_vals, const size_t dimension)
{
	switch (proc_type)
	{
	case eprocess_type::by_rows:
	{
		printf("\nAverage values in rows:\n");
		for (size_t i = 0; i < dimension; ++i)
		{
			printf("Row %u: %lf\n", i, average_vals[i]);
		}
		break;
	}
	case eprocess_type::by_cols:
	{
		printf("\nAverage values in columns:\n");
		for (size_t i = 0; i < dimension; ++i)
		{
			printf("Column %u: %lf\n", i, average_vals[i]);
		}
		break;
	}
	default:
	{
		throw("Incorrect value for parameter 'proc_type' in function PrintAverageVals() call!");
	}
	}
}


int main()
{
	const unsigned ERROR_STATUS = -1;
	const unsigned OK_STATUS = 0;

	unsigned status = OK_STATUS;

	try
	{
		srand((unsigned)time(0));

		const size_t numb_rows = 20;
		const size_t numb_cols = 30;

		double** matrix = new double*[numb_rows];
		for (size_t i = 0; i < numb_rows; ++i)
		{
			matrix[i] = new double[numb_cols];
		}

		double* average_vals_in_rows = new double[numb_rows];
		double* average_vals_in_cols = new double[numb_cols];

		InitMatrix(matrix, numb_rows, numb_cols);

		//PrintMatrix(matrix, numb_rows, numb_cols);

		auto start = std::chrono::system_clock::now();

		std::thread first_thr(FindAverageValues, eprocess_type::by_rows, matrix, numb_rows, numb_cols, average_vals_in_rows);
		std::thread second_thr(FindAverageValues, eprocess_type::by_cols, matrix, numb_rows, numb_cols, average_vals_in_cols);

		first_thr.join();
		second_thr.join();

		auto end = std::chrono::system_clock::now();

		
		//PrintAverageVals(eprocess_type::by_rows, average_vals_in_rows, numb_rows);
		//PrintAverageVals(eprocess_type::by_cols, average_vals_in_cols, numb_cols);

		// Print time
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		std::cout << "Avg time: " << duration.count() / 1000.0 << "ms" << std::endl;


		average_vals_in_rows = new double[numb_rows];
		average_vals_in_cols = new double[numb_cols];

		auto parStart = std::chrono::system_clock::now();
		FindAverageValuesCilk(eprocess_type::by_rows, matrix, numb_rows, numb_cols, average_vals_in_rows);
		FindAverageValuesCilk(eprocess_type::by_cols, matrix, numb_rows, numb_cols, average_vals_in_cols);
		auto parEnd = std::chrono::system_clock::now();


		//PrintAverageVals(eprocess_type::by_rows, average_vals_in_rows, numb_rows);
		//PrintAverageVals(eprocess_type::by_cols, average_vals_in_cols, numb_cols);

		auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(parEnd - parStart);
		std::cout << "Avg time: " << duration1.count() / 1000.0 << "ms" << std::endl;

	}
	catch (std::exception& except)
	{
		printf("Error occured!\n");
		except.what();
		status = ERROR_STATUS;
	}

	return status;
}