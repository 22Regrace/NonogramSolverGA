import tkinter as tk
import sys

CELL_SIZE = 30

def DisplayNonogram(gridHeight, gridWidth, inputStream):
    rows = int(gridHeight)
    cols = int(gridWidth)

    grid = []
    for i in range(rows):
        newRow = []
        for j in range(cols):
            newRow.append(int(inputStream[i * cols + j]))
        grid.append(newRow)
    
    root = tk.Tk()
    root.title("Nonogram Display")

    canvas = tk.Canvas(root, width=cols * CELL_SIZE, height=rows * CELL_SIZE)
    canvas.pack()

    for i, row in enumerate(grid):
        for j, val in enumerate(row):
            if val == 1:
                color = "black" 
            else:
                color = "white"
            
            #top left corner coords for each square
            x0 = j * CELL_SIZE
            y0 = i * CELL_SIZE
            
            #bottom right corner coords for each square
            x1 = x0 + CELL_SIZE
            y1 = y0 + CELL_SIZE
            
            canvas.create_rectangle(x0, y0, x1, y1, fill=color, outline="gray")

    root.mainloop()


def main():
    if len(sys.argv) < 4:
        print("Not enough parameters")
    else:
        DisplayNonogram(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()