#include <stdio.h>
#include <vector>
#include <process.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>

// Comment in if you want to enable Wisdom of Crowds
//#define WOC_ENABLED

// Comment in if you want to save your results to the CSV file
//#define SAVE_RESULTS_TO_CSV

#define GENERATIONS 1000
#define POPULATION_SIZE 250
#define ELITES 15
#define MIN_MUTATION_RATE 0.01
#define MAX_MUTATION_RATE 0.05

#define PRINT_INTERVAL 500 //50

const std::string inputFilePath = "TestData/5x5Puzzle.non";

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

	// Skipping the header
	SkipFileLines(file, 1);

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

// This is a function that looks for blocks of 1s from a row or column (since we need to account for continuity) 
// It returns a vector containing the continuous groups from the row or column so that it can be compared to the actual solution
std::vector<int> GetBlocks(const std::vector<int>& line) {
    std::vector<int> blocks;
    int currentBlock = 0;
    
    for (int i = 0; i < line.size(); i++) {
        if (line[i] == 1) {
            currentBlock++;
        } else { // 0 found
            if (currentBlock > 0) {
                blocks.push_back(currentBlock); // Storing previous block size
                currentBlock = 0;
            }
        }
    }
    
    // Handling case where last block is 1
    if (currentBlock > 0) {
        blocks.push_back(currentBlock);
    }
    
    return blocks;
}

// This is a function that calculates the difference between expected and actual blocks for ROWS (do not need separate functions for 1s and 0s anymore)
double CalculateF1(const std::vector<std::vector<int>>& grid, const std::vector<std::vector<int>>& rowHints) {
    double f1 = 0.0;
    int height = grid.size();
    
    for (int i = 0; i < height; i++) {
        // Get expected blocks from hints
        std::vector<int> expectedBlocks = rowHints[i];
        
        // Get actual blocks from particular solution
        std::vector<int> actualBlocks = GetBlocks(grid[i]);
        
        // Calculate difference between expected and actual blocks
        int maxBlocks = std::max(expectedBlocks.size(), actualBlocks.size());
        for (int j = 0; j < maxBlocks; j++) {
            int expected = 0;
			if (j < expectedBlocks.size()) {
			    expected = expectedBlocks[j];
			}

			int actual = 0;
			if (j < actualBlocks.size()) {
			    actual = actualBlocks[j];
			}
            f1 += abs(expected - actual);
        }
    }
    return f1;
}

// This is a function that calculates the difference between expected and actual blocks for COLUMNS
double CalculateF2(const std::vector<std::vector<int>>& grid, const std::vector<std::vector<int>>& colHints) {
    double f2 = 0.0;
    int height = grid.size();
    int width = grid[0].size();
    
    for (int k = 0; k < width; k++) {
        // Extract the column into a vector
        std::vector<int> column;
        for (int i = 0; i < height; i++) {
            column.push_back(grid[i][k]);
        }
        
        // Get expected blocks from hints
        std::vector<int> expectedBlocks = colHints[k];
        
        // Get actual blocks from current solution
        std::vector<int> actualBlocks = GetBlocks(column);
        
        // Calculate difference
        int maxBlocks = std::max(expectedBlocks.size(), actualBlocks.size());
        for (int j = 0; j < maxBlocks; j++) {
            int expected = 0;
			if (j < expectedBlocks.size()) {
			    expected = expectedBlocks[j];
			}

			int actual = 0;
			if (j < actualBlocks.size()) {
			    actual = actualBlocks[j];
			}
            f2 += abs(expected - actual);
        }
    }
    return f2;
}

// This is a function to calculate the fitness for a given nonogram solution
// It subtracts the total error from the size of the binary matrix
// Hence, a higher fitness indicates a stronger solution
void CalculateFitness(NonogramInstance& nonogram, const NonogramData& data) {
	// Getting height and width
	int N = data.height;
	int M = data.width;

	// Calculate both components
	double f1 = CalculateF1(nonogram.grid, data.rowHints);
	double f2 = CalculateF2(nonogram.grid, data.colHints);

	// Total error
	double f = f1 + f2;

	// Transform to maximization problem
	// Added coefficient of 2 to increase range
	nonogram.fitness = (2 * N * M) - f;
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
void ShowNonogram(const std::vector<std::vector<int>> grid, int gridHeight, int gridWidth, std::string title)
{
	//Create a bitstream from grid to be passed as a parameter
	std::string bitStreamResult = GridToBitStream(grid);

	//Convert it to a parameter that can be used
	std::string parameter = "start /B python ./NonogramGUI.py " + std::to_string(gridHeight) + " "
		+ std::to_string(gridWidth) + " " + bitStreamResult + " " + title;
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

//Prints the best and average fitness of a population of Nonogram Instances
void PrintPopulationStats(const std::vector<NonogramInstance>& pop, int genNumber) {
	float avgFitness = 0.0;
	float bestFitness = 0.0;
	for (const NonogramInstance& instance : pop) {
		avgFitness += instance.fitness;
		if (bestFitness == 0 || bestFitness < instance.fitness) {
			bestFitness = instance.fitness;
		}
	}
	avgFitness /= pop.size();
	std::cout << "Gen " << genNumber + 1 << ": best = " << bestFitness
		<< " avg = " << avgFitness << std::endl;
}

//Crosses over two nonograms by selecting rows from one parent or the other to construct the child
NonogramInstance RowCrossover(const NonogramInstance& parent1, const NonogramInstance& parent2) {
	//Create a child of the same size as the parents
	NonogramInstance child;
	int height = parent1.grid.size();
	int width = parent1.grid[0].size();
	child.grid.resize(height, std::vector<int>(width));

	//Randomly select parent 1 or parent 2 and assign each row to that parent's row
	for (int i = 0; i < height; i++) {
		if (rand() % 2 == 0)
			child.grid[i] = parent1.grid[i];
		else
			child.grid[i] = parent2.grid[i];
	}
	return child;
}

//Mutates an input Nonogram Instance by randomly flipping a value from 0 to 1 or vice versa
void TryFlipMutation(NonogramInstance& instance, double mutationRate) {
	int height = instance.grid.size();
	int width = instance.grid[0].size();

	// Only perform mutation with given probability
	if (((double)rand() / RAND_MAX) < mutationRate) {
		// Pick a random cell
		int i = rand() % height;
		int j = rand() % width;

		//Flips the value in that cell
		instance.grid[i][j] = 1 - instance.grid[i][j];
	}
}

NonogramInstance NonogramSolverGA(const std::vector<NonogramInstance>& initialPop, const NonogramData& data) {

	std::vector<NonogramInstance> pop = initialPop;
	std::vector<NonogramInstance> newPop = {};
	float curMutationRate = MIN_MUTATION_RATE;

	for (int gen = 0; gen < GENERATIONS; gen++) {

		newPop = {};

		//ELITISM
		//Sort the population so that the best individuals can be found easily
		std::vector<NonogramInstance> sortedPop = pop;
		std::sort(sortedPop.begin(), sortedPop.end(), [](const NonogramInstance& a, const NonogramInstance& b) {
			return a.fitness > b.fitness;
			});
		// Push the top n elites into newPop
		for (int i = 0; i < ELITES && i < sortedPop.size(); i++) {
			newPop.push_back(sortedPop[i]);
		}

#ifdef WOC_ENABLED
		//For wisdom of crowds we are going to combine the grid at every point and opt for filled/not depending on majority of experts
		NonogramInstance WOC_ELITE = { {} };
		WOC_ELITE.fitness = 0.0;
		WOC_ELITE.grid.resize(newPop[0].grid.size(), std::vector<int>(newPop[0].grid[0].size(), 0));

		//Go through each of the elites grids 
		for (int i = 0; i < newPop[0].grid.size(); i++) {
			for (int j = 0; j < newPop[0].grid[i].size(); j++) {

				int ExpertsFilledSquare = 0;
				int ExpertsNotFilledSquare = 0;
				for (int k = 0; k < newPop.size(); k++) {

					if (newPop[k].grid[i][j] == 1) {
						ExpertsFilledSquare += 1;
					}
					else {
						ExpertsNotFilledSquare += 1;
					}
				}
				//Fill in the WOC with the majority opinion of the experts
				if (ExpertsFilledSquare >= ExpertsNotFilledSquare) {
					WOC_ELITE.grid[i][j] = 1;
				}
				else {
					WOC_ELITE.grid[i][j] = 0;
				}

			}
		}
		//Only show the GUI for nonogram at the end 
		if ((gen + 1) % GENERATIONS == 0) {
			ShowNonogram(WOC_ELITE.grid, WOC_ELITE.grid.size(), WOC_ELITE.grid[0].size(), "WOC");
		}

		if ((gen + 1) % PRINT_INTERVAL == 0) {
			CalculateFitness(WOC_ELITE, data);
			printf("Fitness of WOC: %f\n", WOC_ELITE.fitness);
		}

		//Add WOC to the ELITES to spawn the next kids
		newPop.push_back(WOC_ELITE);

#endif
		for (int i = ELITES; i < pop.size(); i++) {
			//Select the two parents using a roulette wheel method 
			NonogramInstance parent1 = pop[RouletteWheelSelect(pop)];
			NonogramInstance parent2 = pop[RouletteWheelSelect(pop)];

			//Then create a child by crossing over the parents
			NonogramInstance child = { {} };
			child = RowCrossover(parent1, parent2);

			//Then attempt to apply a mutation to the child and add it to the population
			TryFlipMutation(child, curMutationRate);

			newPop.push_back(child);
		}
		//Set the new population and set each of its instance's fitness
		pop = newPop;
		GetPopFitness(pop, data);

		if ((gen + 1) % PRINT_INTERVAL == 0) {
			PrintPopulationStats(pop, gen);
		}

		//Increase the mutation rate based on the amount of generations passed, or cap it at the max. rate
		curMutationRate = MIN_MUTATION_RATE + (MAX_MUTATION_RATE - MIN_MUTATION_RATE) * (gen / (float)GENERATIONS);

		
	}
	//Get the best individual in the final population to return it
	NonogramInstance bestInstance = { {} };
	for (NonogramInstance& ind : pop) {
		if (bestInstance.fitness == 0 || bestInstance.fitness < ind.fitness)
			bestInstance = ind;
	}
	return bestInstance;
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
	std::vector<NonogramInstance> Population = InitializePopulation(POPULATION_SIZE, nonogramData.height, nonogramData.width);
	GetPopFitness(Population, nonogramData);

	NonogramInstance solution = NonogramSolverGA(Population, nonogramData);

	std::cout << "\nThe fitness of the final solution was " << solution.fitness << ".\n";

	// Printing the real nonogram solution (can get rid of this later)
	printf("\nNonogram Goal State: \n");
	PrintNonogramState(nonogramData.goalState);
	ShowNonogram(nonogramData.goalState, nonogramData.height, nonogramData.width, "Optimal Solution");

#ifdef SAVE_RESULTS_TO_CSV
	resultsFile << "\n\n";
	resultsFile.close();
#endif

	//Needs to be last to execute since anything below will not run until the python window is closed
	ShowNonogram(solution.grid, nonogramData.height, nonogramData.width, "Best Solution");

	return 0;
}

