import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


def produce_graphic(num_r, num_c, grid, dominated, cost,
                    pdf, title, col_del):
#    if title not in ["Failed SubCase"]:#['ContradictionForEverySubCase+ResultsInMinimalityContradiction', 'ContradictionForEverySubCase+HasBroadcastScheme','SubCase']:
#        return


    width = num_c/2
    height = num_r/2 + 1
    fig = plt.figure(figsize=(width, height), constrained_layout=True)

    # define coordinate lists for the dominating set and dominated set
    x_coordinate_dom_set, y_coordinate_dom_set = [], []
    x_coordinate_dominated, y_coordinate_dominated = [], []

    col_index = -1
    for col in range(num_c):
        if col not in col_del:
            col_index += 1
        for row in range(num_r):
            # if vertex is broadcasting
            if grid[row][col] != 0:
                x_coordinate_dom_set.append(col)
                y_coordinate_dom_set.append(row)
                # determine broadcasting neighborhood
                x_diamond = np.array([col, col + grid[row][col],
                                      col, col - grid[row][col], col])
                y_diamond = np.array([row + grid[row][col], row,
                                      row - grid[row][col], row,
                                      row + grid[row][col]])
                # plot neighborhood
                plt.plot(x_diamond, y_diamond, color="r")
            # plot vertices to be dominated if applicable
            if dominated != []:
                # vertex is dominated
                if dominated[row][col] != 0:
                    if col not in col_del:
                        x_coordinate_dominated.append(col_index)
                        y_coordinate_dominated.append(row)

    # plot dominating vertices
    plt.scatter(x_coordinate_dom_set, y_coordinate_dom_set, color="black")
    plt.scatter(x_coordinate_dominated, y_coordinate_dominated, color="blue")

    # draw border
    x_G_corners = np.array([0, num_c - 1, num_c - 1, 0, 0])
    y_G_corners = np.array([0, 0, num_r - 1, num_r - 1, 0])
    plt.plot(x_G_corners, y_G_corners)

    # define ticks and limits of axes
    yticks = [i for i in range(min(y_G_corners) - 1, max(y_G_corners) + 2)]
    xticks = [i for i in range(min(x_G_corners) - 1, max(x_G_corners) + 2)]

    # if induction case
    if col_del != []:
        max_col_id = num_c + len(col_del) + 1
        # only include columns not deleted
        reindex_col = sorted([col for col in range(-1, max_col_id)
                              if col not in col_del])
        plt.xticks(xticks, reindex_col)
    else:
        plt.xticks(xticks)

    plt.yticks(yticks)
    plt.ylim(min(y_G_corners) - 1, max(y_G_corners) + 1)
    plt.xlim(min(x_G_corners) - 1, max(x_G_corners) + 1)
    # correct the yaxis
    plt.gca().invert_yaxis()
    plt.title(title)

    # plot orange rectangle if not dealing with an induction case
    if col_del == []:
        plt.plot([4, num_c-4-1, num_c-4-1, 4, 4],
                 [0, 0, num_r-1, num_r-1, 0], color="orange")
    if col_del != []:
        txt = "Cost: {}, Del Col: {}".format(cost, col_del)
    else:
        txt = "Cost: {}".format(cost)

    plt.xlabel(txt)
    plt.grid(True)
    pdf.savefig(fig)
    plt.close()


# Test
# grid = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#         [0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 2, 0, 0, 0],
#         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0]]
# produce_graphic(3, 14, grid, [], 3, "test.pdf", "thetest")




#pdf.close()

