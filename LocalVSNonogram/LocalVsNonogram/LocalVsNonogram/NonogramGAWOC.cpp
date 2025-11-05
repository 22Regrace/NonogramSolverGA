#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

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

	// Making sure that the data was read in correctly and is not NULL
	if (newNonogram.height == 0 || newNonogram.width == 0 || newNonogram.rowHints.empty() || newNonogram.colHints.empty() || newNonogram.goalState.empty())
	{
		return false;
	}

	return true;
}

int main()
{
	// File poiting to the current nonogram data
	std::string fileName = inputFilePath;

	// Reading in and storing the nonogram data
	NonogramData nonogramData;
	bool inputSuccess = GetFileInfo(nonogramData, fileName);
	if (!inputSuccess)
	{
		printf_s("Something went wrong trying to read file input.\n");
		return -1;
	}


	// Rest of the code goes here


	// Printing the real nonogram solution (can get rid of this later)
	PrintNonogramState(nonogramData.goalState);

	return 0;
}