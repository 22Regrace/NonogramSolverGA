#include <stdio.h>
#include <vector>
#include <process.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>

#define GENERATION_SIZE 5

// Comment in if you want to enable Wisdom of Crowds
//#define WOC_ENABLED

// Comment in if you want to save your results to the CSV file
//#define SAVE_RESULTS_TO_CSV

const std::string inputFilePath = "../../../TestData/db/webpbn/1.non";

/// <summary>
/// The nonogram data collected from the input file
/// </summary>
struct NonogramData
{
	int height = 0;
	int width = 0;
	std::vector<std::vector<int>> rowHints;
	std::vector<std::vector<int>> colHints;
	std::vector<std::vector<int>> goalState;
};

/// <summary>
/// Prints the state of a nonogram (map of 1's and 0's)
/// </summary>
/// <param name="state">State of the nonogram to print</param>
void PrintNonogramState(std::vector<std::vector<int>> state)
{
	for (int i = 0; i < state.size(); i++)
	{
		for (int j = 0; j < state[i].size(); j++)
		{
			printf_s("%d ", state[i][j]);
		}
		printf_s("\n");
	}
}

/// <summary>
/// Skips un-needed lines in the input file
/// </summary>
/// <param name="file">File pointer to move</param>
/// <param name="numLines">How many lines to skip</param>
void SkipFileLines(std::ifstream& file, const int numLines)
{
	std::string buffer;
	for (int i = 0; i < numLines; i++)
	{
		std::getline(file, buffer);
	}
}

/// <summary>
/// Reads in and stores the constraints (hints) from the input file
/// </summary>
/// <param name="hints">Output. The vector to store the row or column hints into</param>
/// <param name="numLines">The number of hints to read</param>
/// <param name="file">The file to read from</param>
void ParseHints(std::vector<std::vector<int>>& hints, const int numLines, std::ifstream& file)
{
	// Skipping the two extra lines before the start of the data
	SkipFileLines(file, 2);

	// Looping through and storing the contraints from the input file
	std::string line;
	for (int i = 0; i < numLines; i++)
	{
		std::getline(file, line);

		std::stringstream ss(line);
		std::string currVal;
		std::vector<int> hint;

		// Pushing each comma seperated value onto a vector
		// Getting the hints for the current row/column
		while (std::getline(ss, currVal, ','))
		{
			hint.push_back(std::stoi(currVal));
		}

		// Pushing the hints for each row/column onto a main vector
		hints.push_back(hint);
	}
}

/// <summary>
/// Reading in and storing the goal state from the input file
/// </summary>
/// <param name="goalState">Output. The correct nonogram solution</param>
/// <param name="file">The input file</param>
void ParseGoal(std::vector<std::vector<int>>& goalState, std::ifstream& file)
{
	// Skipping the empty line before the goal
	SkipFileLines(file, 1);

	std::string line;

	// Getting past the 'goal "' and to the actual data
	size_t start;
	std::getline(file, line);
	start = line.find('"') + 1;
	line.erase(0, start);

	// Storing all of the values in the goal state
	// currIndex holds the index of the current character in the line string
	int currIndex = 0;
	for (int i = 0; i < goalState.size(); i++)
	{
		// Looping through and storing the goal state in a 2D vector matching the nonogram size
		for (int j = 0; j < goalState[i].size(); j++)
		{
			char currChar = line[currIndex++];
			goalState[i][j] = std::atoi(&currChar);
		}
	}
}

/// <summary>
/// Reading in all of the information from the input file
/// </summary>
/// <param name="newNonogram">Output. The nonogramData object being stored</param>
/// <param name="fileName">The input file</param>
/// <returns>Whether or not the data was successfully read in</returns>
bool GetFileInfo(NonogramData& newNonogram, std::string fileName)
{
	std::ifstream file(fileName);
	if (!file.is_open())
	{
		printf_s("Error opening file.\n");
		return false;
	}

	std::string line;

	// Skipping the copyright part of the header
	SkipFileLines(file, 5);

	// Getting the width of the nonogram
	size_t space;
	std::getline(file, line);
	space = line.find(' ');
	newNonogram.width = std::stoi(line.substr(space + 1));

	// Getting the height of the nonogram
	std::getline(file, line);
	space = line.find(' ');
	newNonogram.height = std::stoi(line.substr(space + 1));

	// Getting the constraint data for the rows and columns
	ParseHints(newNonogram.rowHints, newNonogram.height, file);
	ParseHints(newNonogram.colHints, newNonogram.width, file);

	// Getting the goal state data
	newNonogram.goalState.resize(newNonogram.height, std::vector<int>(newNonogram.width, 0));
	ParseGoal(newNonogram.goalState, file);

	file.close();

	// Making sure that the data was read in correctly and is not NULL
	if (newNonogram.height == 0 || newNonogram.width == 0 || newNonogram.rowHints.empty() || newNonogram.colHints.empty() || newNonogram.goalState.empty())
	{
		return false;
	}

	return true;
}

// This is a function to print the header info to the results CSV file before any run data gets added
void InitResultsCSV(std::ofstream& file, const std::string fileName, const int height, const int width)
{
	// Printing the header info for each run of the results
#ifdef WOC_ENABLED
	file << "WOC" << ",";
#else
	file << "Standard GA" << ",";
#endif

	file << height << "x" << width << ",";
	file << fileName << "\n";

	file << "Generation" << ","
		<< "Avg fitness" << ","
		<< "Best fitness" << ","
		<< "Best path"
		<< "\n";
}

// This is a struct to store a particular Nonogram solution
struct NonogramInstance {
	std::vector<std::vector<int>> grid;  // Solution
	double fitness;
};

// --Initialization--

// This is a function to generate a random solution for the population
// i.e. randomly populate a matrix with 0s and 1s
NonogramInstance GenerateSolution(int height, int width) {
	NonogramInstance nonogram;
	nonogram.fitness = 0.0;
	nonogram.grid.resize(height, std::vector<int>(width, 0));

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			nonogram.grid[i][j] = rand() % 2; // Forces remainder to 0 or 1
		}
	}

	return nonogram;
}

// This is a function to initialize a population of random solutions
std::vector<NonogramInstance> InitializePopulation(int size, int height, int width) {
	std::vector<NonogramInstance> population;

	for (int i = 0; i < size; i++) {
		population.push_back(GenerateSolution(height, width));
	}

	return population;
}

// --Fitness Assignment--

// This is a function to calculate the difference between the expected
// number of 1s and actual number of 1s for each *row*
double CalculateF1(const std::vector<std::vector<int>>& grid, const std::vector<std::vector<int>>& rowHints) {
	double f1 = 0.0;
	int height = grid.size();
	int width = grid[0].size();

	for (int i = 0; i < height; i++) {
		// Hint gives expected 1s
		int expectedOnes = 0;
		for (int hint : rowHints[i]) {
			expectedOnes += hint;
		}

		// Count actual 1s from particular solution
		int actualOnes = 0;
		for (int j = 0; j < grid[i].size(); j++) {
			actualOnes += grid[i][j];
		}

		// Get difference
		f1 += abs(expectedOnes - actualOnes);
	}

	return f1;
}

// This is a function to calculate the difference between the expected
// number of 1s and actual number of 1s for each *column*
double CalculateF2(const std::vector<std::vector<int>>& grid, const std::vector<std::vector<int>>& colHints) {
	double f2 = 0.0;
	int height = grid.size();
	int width = grid[0].size();

	for (int i = 0; i < width; i++) {
		// Hint gives expected 1s
		int expectedOnes = 0;
		for (int hint : colHints[i]) {
			expectedOnes += hint;
		}

		// Count actual 1s from particular solution
		int actualOnes = 0;
		for (int j = 0; j < height; j++) {
			actualOnes += grid[j][i];
		}

		// Get difference
		f2 += abs(expectedOnes - actualOnes);
	}

	return f2;
}

// This is a function to calculate the difference between the expected
// number of 0s and actual number of 0s for each *row*
double CalculateF3(const std::vector<std::vector<int>>& grid, const std::vector<std::vector<int>>& rowHints) {
	double f3 = 0.0;
	int height = grid.size();
	int width = grid[0].size();

	for (int i = 0; i < height; i++) {
		// Hint gives expected 1s, so we can deduce 0s
		int expectedOnes = 0;
		for (int hint : rowHints[i]) {
			expectedOnes += hint;
		}
		int expectedZeros = width - expectedOnes;

		// Count actual 1s from particular solution
		int actualOnes = 0;
		for (int j = 0; j < grid[i].size(); j++) {
			actualOnes += grid[i][j];
		}
		int actualZeros = width - actualOnes;

		// Get difference
		f3 += abs(expectedZeros - actualZeros);
	}

	return f3;
}

// This is a function to calculate the difference between the expected
// number of 0s and actual number of 0s for each *column*
double CalculateF4(const std::vector<std::vector<int>>& grid, const std::vector<std::vector<int>>& colHints) {
	double f4 = 0.0;
	int height = grid.size();
	int width = grid[0].size();

	for (int i = 0; i < width; i++) {
		// Hint gives expected 1s, so we can deduce 0s
		int expectedOnes = 0;
		for (int hint : colHints[i]) {
			expectedOnes += hint;
		}
		int expectedZeros = height - expectedOnes;

		// Count actual 1s from particular solution
		int actualOnes = 0;
		for (int j = 0; j < height; j++) {
			actualOnes += grid[j][i];
		}
		int actualZeros = height - actualOnes;

		// Get difference
		f4 += abs(expectedZeros - actualZeros);
	}

	return f4;
}

// This is a function to calculate the fitness for a given nonogram solution
// It subtracts the total error from the size of the binary matrix
// Hence, a higher fitness indicates a stronger solution
void CalculateFitness(NonogramInstance& nonogram, const NonogramData& data) {
	// Getting height and width
	int N = data.height;
	int M = data.width;

	// Calculate all four components
	double f1 = CalculateF1(nonogram.grid, data.rowHints);
	double f2 = CalculateF2(nonogram.grid, data.colHints);
	double f3 = CalculateF3(nonogram.grid, data.rowHints);
	double f4 = CalculateF4(nonogram.grid, data.colHints);

	// Total error
	double f = f1 + f2 + f3 + f4;

	// Transform to maximization problem
	nonogram.fitness = (N * M) - f;
}

// This is a function to calculate the fitness for an entire population of solutions
void GetPopFitness(std::vector<NonogramInstance>& population, const NonogramData& data) {
	for (int i = 0; i < population.size(); i++) {
		CalculateFitness(population[i], data);
	}
}

// --Selection--

// This is a function to implement roulette wheel select
int RouletteWheelSelect(const std::vector<NonogramInstance>& population) {
	// Getting total fitness
	double totalFitness = 0.0;
	for (int i = 0; i < population.size(); i++) {
		totalFitness += population[i].fitness;
	}

	// Generate random value between 0 and totalFitness
	double randomValue = ((double)rand() / RAND_MAX) * totalFitness;

	// Find the corresponding solution
	double cumulativeFitness = 0.0;
	for (int i = 0; i < population.size(); i++) {
		cumulativeFitness += population[i].fitness;
		if (cumulativeFitness >= randomValue) { // Select this individual
			return i;
		}
	}

	return population.size() - 1; // No individual selected
}

/// <summary>
/// This function converts the grid to a bitstream to be used by the python program as a parameter
/// </summary>
/// <param name="state"></param>
/// <returns></returns>
std::string GridToBitStream(std::vector<std::vector<int>> state)
{
	std::string bitStream = "";

	for (int i = 0; i < state.size(); i++)
	{
		for (int j = 0; j < state[i].size(); j++)
		{
			bitStream.append(std::to_string(state[i][j]));

		}
	}
	return bitStream;
}

/// <summary>
/// This function takes in the grid and the height and width of the grid and calls the python GUI to show the nonogram
/// </summary>
/// <param name="grid"></param>
/// <param name="gridHeight"></param>
/// <param name="gridWidth"></param>
void ShowNonogram(const std::vector<std::vector<int>> grid, int gridHeight, int gridWidth)
{
	//Create a bitstream from grid to be passed as a parameter
	std::string bitStreamResult = GridToBitStream(grid);

	//Convert it to a parameter that can be used
	std::string parameter = "python ./NonogramGUI.py " + std::to_string(gridHeight) + " "
		+ std::to_string(gridWidth) + " " + bitStreamResult;
	int result = system(parameter.c_str());
}

// This function updates the results CSV file to print:
// generation #, average fitness, bestfitness, and best path in the generation
void UpdateResultsCSV(std::ofstream& file, int gen, double avgFitness, NonogramInstance bestGrid)
{
	file << gen << ","
		<< avgFitness << ","
		<< bestGrid.fitness << ",";

	for (int i = 0; i < bestGrid.grid.size(); i++)
	{
		for (int j = 0; j < bestGrid.grid[i].size(); j++)
		{
			file << bestGrid.grid[i][j];
		}
	}
	file << "\n";
}

int main()
{
	// File poiting to the current nonogram data
	std::string fileName = inputFilePath;

	// Generating random seed for solution generation
	// Process-specific value ensures uniqueness
	unsigned long long seed = (unsigned long long)time(0) ^ _getpid();
	srand(seed);
	
	// Reading in and storing the nonogram data
	NonogramData nonogramData;
	bool inputSuccess = GetFileInfo(nonogramData, fileName);
	if (!inputSuccess)
	{
		printf_s("Something went wrong trying to read file input.\n");
		return -1;
	}

#ifdef SAVE_RESULTS_TO_CSV
	// Opening the results CSV file
	std::ofstream resultsFile("NonogramSolverResults.csv", std::ios::out | std::ios::app);
	if (!resultsFile.is_open())
	{
		printf_s("Error opening file to write results.\n");
		return -1;
	}
	InitResultsCSV(resultsFile, fileName, nonogramData.height, nonogramData.width);
#endif

	// Rest of the code goes here
	std::vector<NonogramInstance> Population = InitializePopulation(GENERATION_SIZE, nonogramData.height, nonogramData.width);
	GetPopFitness(Population, nonogramData);

	for (int i = 0; i < GENERATION_SIZE; i++) 
	{
		ShowNonogram(Population[i].grid, nonogramData.height, nonogramData.width);
	}

	// Printing the real nonogram solution (can get rid of this later)
	printf("\nSolution Nonogram: \n");
	PrintNonogramState(nonogramData.goalState);

#ifdef SAVE_RESULTS_TO_CSV
	resultsFile << "\n\n";
	resultsFile.close();
#endif

	return 0;
}

