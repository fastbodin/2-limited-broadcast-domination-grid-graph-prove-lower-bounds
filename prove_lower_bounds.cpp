#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gurobi_c++.h"

////****** SET PROBLEM INPUTS ******//
//// number of rows in the graph
#define NUM_R 6
//// number of columns in the graph
//// NOTE: NUM_C >= 13
#define NUM_C 13
//// set to 1 if the graph is cycle x cycle
//// set to 0 if the graph is path x cycle
#define CYCLE 1
//// minimum cost of a broadcast
#define MIN_C 2
//// maximum cost of a broadcast
#define MAX_C 3
//// do you want to produce figures?
#define PROD_FIG 1
//// induction values
//// first index needs to be zero
//// for C_3 x C_n
//// int mvalues[] = {0, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 8, 8, 9, 10};
int mvalues[] = {0, 2, 3, 4, 5, 7, 7, 9, 10, 11, 12, 14, 14, 16, 17, 18, 18};
////****** SET PROBLEM INPUTS ******//

//****** COUNTER ******//
//
// for each possible cost of broadcast, keep a counter to 
// track things
struct Counter{
    // number of cases at this cost considered
    unsigned long long int C;
    // those which do not dominate the middle
    unsigned long long int DoesNotDominate;
    // those which contain forbidden broadcasts
    unsigned long long int ForbiddenBroadcast; 
    // those which have base broadcast contradictions
    unsigned long long int HasBroadcast;
    // those which result in a minimality contradiction
    unsigned long long int InductiveArgument;
    // those which result in a neccessary broadcast contradiction 
    // + HasBroadcast contradiction
    unsigned long long int NecessaryBroadcastHasBroadcast;
    // those which result in a neccessary broadcast contradiction 
    // + InductiveArgument contradiction 
    unsigned long long int NecessaryBroadcastInductiveArgument;
    // those which have a subcase contradiction
    // + HasBroadcast contradiction
    unsigned long long int AllSubcasesHasBroadcast;
    // those which have a subcase contradiction
    // + InductiveArgument contradiction
    unsigned long long int AllSubcasesInductiveArgument; 
    // case which failed
    int Fail;
};
//
//****** COUNTER ******//


//****** PRODUCE FIGURES ******//
//
// this is for producing figures
void ProduceGraphic() {
    printf("produce_graphic(%d, num_c, grid, dominated, cost, pdf, title, col_del)\n", NUM_R);
};

// produce the broadcast information for the figures
void ProduceGrid(int grid[NUM_R][NUM_C], int cost) {
    int col, row;
    printf("grid = [[");
    for (row = 0; row < NUM_R; row++) {
        for (col = 0; col < NUM_C; col++) {
            if (col == NUM_C -1) {
                if (row == NUM_R-1) {
                    printf("%d]]\n", grid[row][col]);
                } else {
                    printf("%d], [", grid[row][col]);
                }
            } else {
                printf("%d, ", grid[row][col]);
            }
        }
    }
    printf("cost = %d\n", cost);
};

// produce the dominating set for the figures
void ProduceDominated(int dominated[NUM_R][NUM_C]) {
    int col, row;
    printf("dominated = [[");
    for (row = 0; row < NUM_R; row++) {
        for (col = 0; col < NUM_C; col++) {
            if (col == NUM_C-1) {
                if (row == NUM_R-1) {
                    printf("%d]]\n", dominated[row][col]);
                } else {
                    printf("%d], [", dominated[row][col]);
                }
            } else {
                printf("%d, ", dominated[row][col]);
            }
        }
    }
};

// produce the columns which are deleted for the fgures
void ProduceColDel(int num_del_col, int col_del_order[NUM_C]) {
    int col, row;
    printf("col_del = [");
    for (col = 0; col < num_del_col; col++) {
        if (col == num_del_col -1){
            printf("%d]\n", col_del_order[col]);
        } else {
            printf("%d, ", col_del_order[col]);
        }
    }
};
//
//****** PRODUCE FIGURES ******//


//****** PRODUCE VERTEX BROADCAST INFORMATION ******//
//
// we define a struct for each vertex in the graph
// this struct contains the information about which 
// vertices it dominates at given strengths and distances
struct Vertex{
    // dist[str] contains a list of the vertices (stored by [row, col])
    // which are dominated by the given vertex when broadcasting at 
    // strength (str + 1)
    // e.g. if dist[1][0][0], dist[1][0][1] = 2, 3 then the vertex in row
    // 2 and column 3 is dominated if the given vertex broadcasts at 
    // strength 2
    //
    // len[str] = the number of vertices dominated by the given vertex
    // if broadcasting at strength (str + 1)
    //
    // note that a vertex broadcasting at strength 2 (and therefore 1)
    // can dominate at most 13 other vertices in path x cycle and
    // cycle x cycle
    int dist[2][13][2], len[2];
    
    // dist_E[str] contains a list of the vertices (stored by [row, col])
    // which are at EXACTLY distance (str + 1) from the given vertex and 
    // are therefore dominated by the given vertex if broadcasting 
    // at strength (str + 1)
    //
    // len_E[str + 1] = the number of vertices dominated by the given
    // vertex (at distance EXACTLY (str + 1)) if broadcasting at 
    // strength (str + 1)
    // 
    // note that the number of vertices at distance EXACTLY 2
    // from a given vertex is 8
    // note that the number of vertices at distance EXACTLY 1
    // from a given vertex is 4
    int dist_E[2][8][2], len_E[2];

    // the following information is used when considering possible
    // subcases. when creating all possible subcases, we iteratively
    // add broadcast vertices in columns < 4 and > NUM_C - 5 which
    // dominate vertices in columns 4 and NUM_C - 5 (resp.)

    // if a vertex is in column 4 or NUM_C-5, we keep a list
    // of all the vertices in columns < 4 and > NUM_C-5 (resp.)
    // which dominate it when broadcasting at different strengths
    // e.g. if dom_by[str][1][0], dom_by[str][1][1] = 2, 3 then
    // the given vertex is dominated by the vertex in row 2 and 
    // columns 3 when broadcasting at strength (str + 1)
    //
    // num_dom_by[str] = the number of vertices in columns < 4
    // or > NUM_C-5 (dependant upon the given vertex) which dominate
    // it when broadcasting at strength (str + 1)
    //
    // note it can be dominated by at most 4 vertices from such columns
    int dom_by[2][4][2], num_dom_by[2];

    // if a vertex is in columns < 4 or > NUM_C-5, we keep a
    // bit packed integer of the vertices it dominates in columns 4 and 
    // NUM_C-5 (resp.) when broadcasting at strength 1 and 2
    // e.g. if the given vertex dominates a vertex in row 3 of columns
    // NUM_C-5 when broadcasting at strength 2 then the third bit of
    // dom_col_NUMC_5[1] is a 1.
    uint32_t dom_col_4[2];
    uint32_t dom_col_NUMC_5[2];
};
//
//****** PRODUCE VERTEX BROADCAST INFORMATION ******//


//****** SOME BASIC FUNCTIONS WHICH ARE HELPFUL ******//
//
// output the minimum of two integers
int min(int val_1, int val_2) {
    if (val_1 <= val_2) {
        return val_1;
    }
    return val_2;
};

// determine the distance between two vertices 
int distance(int vrow, int vcol, int urow, int ucol) {
    // if we are dealing with a cycle x cycle
    if (CYCLE) {
        return abs(vcol - ucol) + min(abs(vrow - urow), NUM_R - abs(vrow - urow));
    }
    // else path x cycle
    return abs(vcol - ucol) + abs(vrow - urow);
}
 
// update broadcast (stored in grid) to include the vertex in (row, col) 
// broadcasting at strength str
void add_broadcast(int row, int col, int str, int grid[NUM_R][NUM_C],
                   int dominated[NUM_R][NUM_C],
                   int num_v, int v_dom[13][2]) {
    int i;

    // vertex is now broadcasting at strength str
    grid[row][col] = str;
    // for each vertex is dominates with given strength str
    for (i = 0; i < num_v; i++) {
        // update dominated array
        // v_dom[i][0] = row of vertex i
        // v_dom[i][1] = col of vertex i
        dominated[v_dom[i][0]][v_dom[i][1]] += str;
    }
};

// update broadcast (stored in grid) to include the vertex in (row, col) 
// broadcasting at strength 0
void rem_broadcast(int row, int col, int str, int grid[NUM_R][NUM_C],
              int dominated[NUM_R][NUM_C], int num_v, int v_dom[13][2]) {
    int i;

    // vertex is now broadcasting at strength 0
    grid[row][col] = 0;
    // for each vertex is dominates with given strength str
    for (i = 0; i < num_v; i++) {
        // update dominated array
        // v_dom[i][0] = row of vertex i
        // v_dom[i][1] = col of vertex i
        dominated[v_dom[i][0]][v_dom[i][1]] -= str;
    }
};
//
//****** SOME BASIC FUNCTIONS WHICH ARE HELPFUL ******//

   
//****** FILL IN THE INFORMATION IN THE VERTEX STRUCTS ******//
//
// graph[i][j][k] contains the Vertex struct of the vertex 
// in row j and column k of the graph after having deleted i
// columns from the graph in the induction step
void make_graph(struct Vertex graph[NUM_C][NUM_R][NUM_C]) {
    int vrow, vcol, urow, ucol, dist, num_c, num_del_col, str, v_id;
    // pointer for vertex
    struct Vertex *ptr;

    // THE FOLLOWING LOOP FILLS IN THE INFORMATION FOR DIST, LEN,
    // DIST_E, AND LEN_E FOR EACH VERTEX
    //
    // for each possible number of deleted columns
    for (num_del_col = 0; num_del_col < NUM_C; num_del_col++){
        // after having deleted num_del_col colums, the graph
        // has NUM_C - num_del_col columns left
        // note: offset by one to get the correct index
        num_c = NUM_C - 1 - num_del_col;
        // for each possible strength (1 and 2)
        for (str = 1; str < 3; str++){
            // iterate over all the vertices
            for (vrow = 0; vrow < NUM_R; vrow++){
                for (vcol = 0; vcol <= num_c; vcol++){                
                    // add a pointer to vertex
                    ptr = &graph[num_del_col][vrow][vcol];
                    // iterate over all possible vertices it might
                    // dominate
                    for (urow = 0; urow < NUM_R; urow++){
                        for (ucol = 0; ucol <= num_c; ucol++){
                            // determine distance between the two vertices
                            dist = distance(vrow, vcol, urow, ucol);
                            // if distance is within strength
                            if (dist <= str) {
                                // add information that vertex dominates
                                // the vertex in row urow and col ucol
                                ptr->dist[str - 1][ptr->len[str - 1]][0] = urow;
                                ptr->dist[str - 1][ptr->len[str - 1]][1] = ucol;
                                // incrememnt the number of vertices it
                                // is known to dominate
                                ptr->len[str - 1] += 1;
                            }
                            // if distance is exactly strength
                            if (dist == str){
                                // add information that vertex dominates
                                // the vertex in row urow and col ucol
                                ptr->dist_E[str - 1][ptr->len_E[str - 1]][0] = urow;
                                ptr->dist_E[str - 1][ptr->len_E[str - 1]][1] = ucol;
                                // incrememnt the number of vertices it
                                // is known to dominate
                                ptr->len_E[str - 1] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    // THE FOLLOWING LOOP FILLS IN THE INFORMATION FOR DOM_BY, NUM_DOM_BY,
    // DOM_COL_4, AND DOM_COL_NUMC_5 IF THE VERTEX IS IN THE CORRECT
    // COLUMNS
    //
    // for each possible strength (1 and 2)
    for (str = 1; str < 3; str++) {
        // for the vertices in column 4
        for (vrow = 0; vrow < NUM_R; vrow++) {
            // iterate over vertices in columns < 4
            for (urow = 0; urow < NUM_R; urow++) {
                for (ucol = 0; ucol < 4; ucol++) {
                     // determine distance between the two vertices
                    dist = distance(vrow, 4, urow, ucol);                   
                    if (dist <= str) {
                        // add information that vertex in row vrow and
                        // column 4 is dominated by
                        // the vertex in row urow and column ucol
                        ptr = &graph[0][vrow][4];
                        ptr->dom_by[str - 1][ptr->num_dom_by[str - 1]][0] = urow;
                        ptr->dom_by[str - 1][ptr->num_dom_by[str - 1]][1] = ucol;
                        // incrememnt the number of vertices it
                        // is known to be dominated by in col < 4
                        ptr->num_dom_by[str - 1] += 1;
                        // the bit corresponding to vrow is changed to a 1
                        // since the vertex in row urow and column ucol
                        // dominates it when broadcasting at strength (str - 1)
                        //
                        // note that we only keep this information for the 
                        // graph after having deleted no columns
                        graph[0][urow][ucol].dom_col_4[str - 1] |= 1 << vrow;

                    }
                }
                // iterate over vertices in columns > NUM_C-5
                for (ucol = NUM_C-4; ucol < NUM_C; ucol++) {
                     // determine distance between the two vertices
                    dist = distance(vrow, NUM_C-5, urow, ucol);                   
                    if (dist <= str) {
                        // add information that vertex in row vrow and
                        // column NUM_C-5 is dominated by
                        // the vertex in row urow and column ucol
                        ptr = &graph[0][vrow][NUM_C-5];
                        ptr->dom_by[str - 1][ptr->num_dom_by[str - 1]][0] = urow;
                        ptr->dom_by[str - 1][ptr->num_dom_by[str - 1]][1] = ucol;
                        // incrememnt the number of vertices it
                        // is known to be dominated by in col > NUM_C-5
                        ptr->num_dom_by[str - 1] += 1;
                        // the bit corresponding to vrow is changed to a 1
                        // since the vertex in row urow and column ucol
                        // dominates it when broadcasting at strength (str - 1)
                        //
                        // note that we only keep this information for the 
                        // graph after having deleted no columns
                        graph[0][urow][ucol].dom_col_NUMC_5[str - 1] |= 1 << vrow;
                    }
                }
            }
        }
    }
};
//
//****** FILL IN THE INFORMATION IN THE VERTEX STRUCTS ******//


//****** LP TOOLS ******//
//
// after having deleted num_del_col columns, is there a broadcast
// which dominates the remaining vertices of cost <= x?
bool HasBroadcast(int num_del_col, int col_del_order[NUM_C],
                        int dominated[NUM_R][NUM_C],
                        struct Vertex graph[NUM_R][NUM_C],
                        int x, GRBEnv grb_env) {
    int row, col, v, num_c, i, col_not_del;
    int col_index = 0;
    // pointer for vertex information
    struct Vertex *ptr;

    // create an empty model
    GRBModel model = GRBModel(grb_env);
    // define a 2-dim array for the variables
    // vars[i][row][col] = variables for vertex in [row][col] broadcasting
    // strength i + 1
    // define number of columns and number of variables
    num_c = NUM_C - num_del_col;
    GRBVar vars[2][NUM_R][num_c];

    // make variables
    for (col = 0; col < num_c; col++) {
        for (row = 0; row < NUM_R; row++) {
            // lower bound, upper bound, objective coeff, type
            vars[0][row][col] = model.addVar(0, 1, 1, GRB_BINARY);
            vars[1][row][col] = model.addVar(0, 1, 2, GRB_BINARY);
        }
    }

    // if no column has been deleted from the graph
    if (num_del_col == 0) {
        // define constraints
        for (col = 0; col < num_c; col++) {
            for (row = 0; row < NUM_R; row++) {
                // if the vertex is dominated in the broadcast we hope to
                // beat, it must be dominated in the new broadcast
                if (dominated[row][col] != 0) {
                    GRBLinExpr expr = 0;
                    ptr = &graph[row][col];
                    // consider vertices within distance 1
                    for (v = 0; v < ptr->len[0]; v++) {
                        // add the variables corresonding to said vertex
                        // broadcasting at strength 1
                        expr += vars[0][ptr->dist[0][v][0]][ptr->dist[0][v][1]];
                    }
                    // consider vertices within distance 2
                    for (v = 0; v < ptr->len[1]; v++) {
                        // add the variables corresonding to said vertex
                        // broadcasting at strength 2
                        expr += vars[1][ptr->dist[1][v][0]][ptr->dist[1][v][1]];
                    }
                    // add the constraint to the model
                    model.addConstr(expr >= 1);
                }
            }
        }
    // columns have been deleted from the graph
    } else {
        // determine the indices of the columns which have not been deleted
        for (col = 0; col < NUM_C; col++) {
            // first assume the column has not been deleted
            col_not_del = 1;
            // check over the list of deleted columns
            for (i = 0; i < num_del_col; i++) {
                // if we want to delete the columns
                if (col == col_del_order[i]) {
                    col_not_del = 0;
                    break;
                }
            }
            // for the columns that have not been deleted
            if (col_not_del) {
                // for the vertices in this columns
                for (row = 0; row < NUM_R; row++) {
                    // that are dominated in the previous broadcast
                    if (dominated[row][col] != 0) {
                        GRBLinExpr expr = 0;
                        // col_index = the index of the next column
                        // we will consider when adding constraints
                        // e.g. if column 2 was deleted and columns 0, 1,
                        // and 3 were not, then (from the variables 
                        // persective) column 3 is really just column 2
                        // i.e., we stitch the graph back together
                        // after having deleted column 2
                        ptr = &graph[row][col_index];
                        // consider vertices within distance 1
                        for (v = 0; v < ptr->len[0]; v++) {
                            // add the variables corresonding to said vertex
                            // broadcasting at strength 1
                            expr += vars[0][ptr->dist[0][v][0]][ptr->dist[0][v][1]];
                        }
                        // consider vertices within distance 2
                        for (v = 0; v < ptr->len[1]; v++) {
                            // add the variables corresonding to said vertex
                            // broadcasting at strength 2
                            expr += vars[1][ptr->dist[1][v][0]][ptr->dist[1][v][1]];
                        }
                        // add the constraint to the model
                        model.addConstr(expr >= 1);
                    }
                }
                // column was added
                col_index += 1;
            }
        }
    }

    // The objective is to minimize the costs
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // The objective coefficients are set during the creation of
    // the decision variables above. Run model
    model.optimize();

    // does there exist a broadcast of cost <= x?
    if (int(model.get(GRB_DoubleAttr_ObjVal)) <= x) {
        //****** FIGURES ******//
        if (PROD_FIG) {
            // determine what the contradictory broadcast looks like
            int new_grid[NUM_R][NUM_C] = {0};
            int new_cost = int(model.get(GRB_DoubleAttr_ObjVal));
            for (row = 0; row < NUM_R; row++) {
                for (col = 0; col < num_c; col++) {
                    new_grid[row][col] = int(vars[0][row][col].get(GRB_DoubleAttr_X)) +
                                         2*int(vars[1][row][col].get(GRB_DoubleAttr_X));
                }
            }        
            ProduceGrid(new_grid, new_cost);
        }
        //****** FIGURES ******//
        //
        // found a contradiction
        return 1;
    }
    // did not find a contradiction
    return 0;
};

// does the current broadcast result in a minimality contradiction?
bool InductiveArgument(int dominated[NUM_R][NUM_C],
                                      struct Vertex graph[NUM_C][NUM_R][NUM_C],
                                      int x, GRBEnv grb_env) {
    int col_del_index = 0;
    int col_del_order[NUM_C];
    int num_dom_per_col[NUM_C] = {0};
    int col, row, num_v, num_del_col;

    // determine how many vertices are dominated in each col
    for (col = 0; col < NUM_C; col++) {
        for (row = 0; row < NUM_R; row++) {
            // if vertex is dominated
            if (dominated[row][col] != 0) {
                num_dom_per_col[col] += 1;
            }
        }
    }
    // based on how many vertices they dominate, make 
    // an ordered list of the columns which dominate a 
    // non-zero number of vertices
    for (num_v = NUM_C; num_v > 0; num_v--) {
        for (col = 0; col < NUM_C; col++) {
            // if column has num_v vertices dominated
            if (num_dom_per_col[col] == num_v) {
                col_del_order[col_del_index] = col;
                col_del_index += 1;
            }
        }
    }
    // try deleting columns for induction contradiction
    // there are col_del_index + 1 columns with a non-zero number
    // of dominated vertices, only need to consider them
    for (num_del_col = 1; num_del_col < col_del_index + 1; num_del_col++) {
        // x - mvalues[num_del_col] = cost we want to beat given
        // by the induction values
        // is there such a contradictory broadcast?
        if (HasBroadcast(num_del_col, col_del_order, dominated,
                                graph[num_del_col], x - mvalues[num_del_col],
                                grb_env)) {
            //****** FIGURES ******//
            if (PROD_FIG) {
                ProduceColDel(num_del_col, col_del_order);
                printf("num_c = %d\n", NUM_C - num_del_col);
            }
            //****** FIGURES ******//
            //
            // found a contradiction
            return 1;
        }
    }
    // tried deleting all the possible number of columns and never
    // found a contradiction
    return 0;
};
//
//****** LP TOOLS ******//


//****** USEFUL TOOLS TO HAVE WHEN CREATING THE NECCESSARY BROADCAST ******//
//
// add neccessary broadcast vertices
int AddNecBroadcast(int grid[NUM_R][NUM_C], int dominated[NUM_R][NUM_C],
                    struct Vertex graph[NUM_R][NUM_C], int cost) {
    int row;
    int nec_cost = cost;

    for (row = 0; row < NUM_R; row++) {
        // if there is a vertex in column 5 not dominated
        if (dominated[row][5] == 0) {
            // add vertex in same row, column 3 broadcasting at strength 2
            add_broadcast(row, 3, 2, grid, dominated, graph[row][3].len[1], 
                            graph[row][3].dist[1]);
            // the cost goes up
            nec_cost += 2;
        }
        // if there is a vertex in column NUM_C-6 not dominated
        if (dominated[row][NUM_C-6] == 0) {
            // add vertex in same row, column NUM_C-4 broadcasting at strength 2
            add_broadcast(row, NUM_C-4, 2, grid, dominated, graph[row][NUM_C-4].len[1], 
                            graph[row][NUM_C-4].dist[1]);
            // cost goes up
            nec_cost += 2;
        }
    }
    // return the new cost of the neccessary broadcast
    return nec_cost;
};

// remove neccessary broadcast vertices
void RemoveNecBroadcast(int grid[NUM_R][NUM_C], int dominated[NUM_R][NUM_C],
                            struct Vertex graph[NUM_R][NUM_C]) {
    int row;
    for (row = 0; row < NUM_R; row++) {
        if (grid[row][3] == 2) {
            // add vertex in same row, column 3 broadcasting at strength 2
            rem_broadcast(row, 3, 2, grid, dominated, graph[row][3].len[1], 
                            graph[row][3].dist[1]);
        }
        if (grid[row][NUM_C-4] == 2) {
            // add vertex in same row, column NUM_C-4 broadcasting at strength 2
            rem_broadcast(row, NUM_C-4, 2, grid, dominated, graph[row][NUM_C-4].len[1], 
                            graph[row][NUM_C-4].dist[1]);
        }
    }
};
//
//****** USEFUL TOOLS TO HAVE WHEN CREATING THE NECCESSARY BROADCAST ******//


//****** USEFUL TOOLS TO HAVE WHEN CREATING THE SUBBROADCASTS ******//
//
// check if there would be a forbidden broadcast 
// when adding a vertex of strength 1
bool NoLocalForbiddenBroadcastStrength1(int grid[NUM_R][NUM_C],
                             struct Vertex vertex) {
    int i;
    // if it has a neighbor broadcasting at non-zero strength
    for (i = 0; i < vertex.len_E[0]; i++) {
        if (grid[vertex.dist_E[0][i][0]][vertex.dist_E[0][i][1]]!= 0) {
            return 0;
        }
    }
    // if it has a vertex at distance 2 broadcasting at strength 1
    for (i = 0; i < vertex.len_E[1]; i++) {
        if (grid[vertex.dist_E[1][i][0]][vertex.dist_E[1][i][1]] == 1) {
            return 0;
        }
    }
    // no forbidden broadcast found
    return 1;
};

// check if there would be a forbidden broadcast
// when adding a vertex of strength 2
bool NoLocalForbiddenBroadcastStrength2(int grid[NUM_R][NUM_C],
                             struct Vertex vertex) {
    int i;
    // if it has a neighbor broadcasting at strength 1
    for (i = 0; i < vertex.len_E[0]; i++) {
        if (grid[vertex.dist_E[0][i][0]][vertex.dist_E[0][i][1]] == 1) {
            return 0;
        }
    }
    // no forbidden broadcast found
    return 1;
};

// if we want to add a new vertex broadcasting in the subcase, does it
// dominate (in column 4 or NUM_C-5, resp.) a superset 
// (of the vertices in column 4 or NUM_C-5, resp.) of 
// some other vertex which we previously added to the broadcast?
bool NoDominationOverlap(int vcol, struct Vertex new_vertex, int str,
                            int num_added[2], uint32_t dom_sets[2][NUM_R]) {
    // num_added[0] = number of vertices which are added to
    // the broadcast to dominate column 4
    // num_added[1] = number of vertices which are added to
    // the broadcast to dominate column NUM_C-5
    //
    // dom_sets[0][i] = bit-packed integers corresonding to the rows of the
    // vertices dominated in column 4 by the ith vertex added to dominate column 4
    // dom_sets[1][i] = bit-packed integers corresonding to the rows of the
    // vertices dominated in column NUM_C-5 by the ith vertex added to
    // dominate column NUM_C-5
    int i;
    // if we are dealing with vertices which dominate column 4
    if (vcol < 4) {
        // if there are no vertices already dominating column 4
        if (num_added[0] == 0) {
            // cannot be a superset then
            return 1;
        }
        // check if the new dominating set is a superset of previously 
        // added dominating sets
        for (i = 0; i < num_added[0]; i++) {
            if ((new_vertex.dom_col_4[str - 1] & dom_sets[0][i]) == dom_sets[0][i]) {
                // if it is, we wont add this vertex
                return 0;
            }
        }
        // no old dominating set is a subset of the new dominating set (when restricted
        // to column 4)
        return 1;
    } else {
        // if there are no vertices already dominating column NUM_C-5
        if (num_added[1] == 0) {
            // cannot be a superset then
            return 1;
        }
        // check if the new dominating set is a superset of previously 
        // added dominating sets
        for (i = 0; i < num_added[1]; i++) {
            if ((new_vertex.dom_col_NUMC_5[str - 1] & dom_sets[1][i]) == dom_sets[1][i]) {
                // if it is, we wont add this vertex
                return 0;
            }
        }
        // no old dominating set is a subset of the new dominating set (when restricted
        // to column NUM_C-5)
        return 1;
    }
};

// update the number of vertices added to dominate vertices in column 4 or
// NUM_C-5 and also update the list which holds the bit-packed integers
// which indicate which vertices (by row) they dominate
void add_to_dom_set_and_num_added(int vcol, struct Vertex v_ptr, int str,
                                int num_added[2], uint32_t dom_sets[2][NUM_R]) {
    if (vcol < 4) {
        dom_sets[0][num_added[0]] = v_ptr.dom_col_4[str - 1];
        num_added[0] += 1;
    } else {
        dom_sets[1][num_added[1]] = v_ptr.dom_col_NUMC_5[str - 1];
        num_added[1] += 1;
    }
}

// update the number of vertices added to dominate vertices in column 4
// or NUM_C-5 since we removed a vertex
void rem_from_dom_set_and_num_added(int vcol, int num_added[2]) {
    if (vcol < 4) {
        num_added[0] -= 1;
    } else {
        num_added[1] -= 1;
    }
}
//
//****** USEFUL TOOLS TO HAVE WHEN CREATING THE SUBBROADCASTS ******//


//****** CONSIDER ALL POSSIBLE SUBBROADCASTS ******//
//
// recursion to produce all possible subbroadcasts
bool ProduceAllSubBroadcasts(int level, int grid[NUM_R][NUM_C], int dominated[NUM_R][NUM_C],
                             struct Vertex graph[NUM_C][NUM_R][NUM_C], int cost, 
                             struct Counter counter[MAX_C + 1], int undom_verts[2*NUM_R][2],
                             int num_undom, GRBEnv grb_env, int num_added[2],
                             uint32_t dom_sets[2][NUM_R], int counter_cost) {

    int row, col, v_id, vrow, vcol;
    struct Vertex *ptr;
    struct Vertex *v_ptr;

    // if we were unable to prove a case
    if (counter[counter_cost].Fail != 0) {
        return 0;
    }

    // if the solution is now unique (i.e. all the vertices that remained to be
    // dominated are dominated)
    if (level >= num_undom) {
        //****** FIGURES ******//
        if (PROD_FIG) {
            printf("dominated = []\n");
            ProduceGrid(grid, cost);
            printf("title = 'SubCase'\n");
            printf("col_del = []\n");
            printf("num_c = %d\n", NUM_C);
            ProduceGraphic();
            ProduceDominated(dominated);
        }
        //****** FIGURES ******//
        //
        // first we test whether or not the vertices dominated by the 
        // assumed subbroadcast can be dominated with less cost
        if (HasBroadcast(0, NULL, dominated, graph[0], cost - 1, grb_env)) {
            //****** FIGURES ******//
            if (PROD_FIG) {
                printf("col_del = []\n");
                printf("num_c = %d\n", NUM_C);
                printf("title = 'ContradictionForEverySubCase+HasBroadcast'\n");
                ProduceGraphic();
            }
            //****** FIGURES ******//
            //
            // we found a contradiction
            counter[counter_cost].AllSubcasesHasBroadcast += 1;
            return 1;
        }
        // the vertices dominated by the assumed subbroadcast does not yield a
        // contradiction as is, try induciton
        if (InductiveArgument(dominated, graph, cost, grb_env)) {
            //****** FIGURES ******//
            if (PROD_FIG) {
                printf("title = 'ContradictionForEverySubCase+InductiveArgument'\n");
                ProduceGraphic();
            }
            //****** FIGURES ******//
            //
            // we found a conradiction
            counter[counter_cost].AllSubcasesInductiveArgument += 1;
            return 1;
        }
        //****** FIGURES ******//
        //if (PROD_FIG) {
        ProduceGrid(grid, cost);
        printf("col_del = []\n");
        printf("num_c = %d\n", NUM_C);
        printf("dominated = []\n");
        printf("title = 'Failed SubCase'\n");
        ProduceGraphic();
        //}
        //****** FIGURES ******//
        //
        // could not find a contradiction for this subcase
        counter[0].Fail += 1;
        counter[counter_cost].Fail += 1;
        return 0;
    }

    // define the row and column of next vertex we want to dominate (if
    // not already dominated)
    row = undom_verts[level][0], col = undom_verts[level][1];
    // is the vertex dominated?
    if (dominated[row][col] != 0) {
        ProduceAllSubBroadcasts(level + 1, grid, dominated, graph, cost,
                                counter, undom_verts, num_undom, grb_env,
                                num_added, dom_sets, counter_cost);
        // no more options
        return 0;
    }
    // vertex is not dominated
    ptr = &graph[0][row][col];
    // try adding a broadcasting at strength 2
    // note we consider [1] since we are adding strength 2
    for (v_id = 0; v_id < ptr->num_dom_by[1]; v_id++) {
        // define coordinates of the vertex we want to add 
        // which dominates (row, col)
        vrow = ptr->dom_by[1][v_id][0], vcol = ptr->dom_by[1][v_id][1];
        // if vertex is not already broadcasting
        if (grid[vrow][vcol] == 0) {
            v_ptr = &graph[0][vrow][vcol];
            // if vertex does not create a forbidden broadcast
            if (NoLocalForbiddenBroadcastStrength2(grid, *v_ptr)) { 
                // if there is not overlap in the broadcast by this vertex
                if (NoDominationOverlap(vcol, *v_ptr, 2, num_added, dom_sets)) {
                    // add this vertex to the broadcast
                    add_broadcast(vrow, vcol, 2, grid, dominated, v_ptr->len[1], v_ptr->dist[1]);
                    // update the num_added and dom_sets
                    add_to_dom_set_and_num_added(vcol, *v_ptr, 2, num_added, dom_sets);
                    // next step of recursion
                    ProduceAllSubBroadcasts(level + 1, grid, dominated, graph,
                                            cost + 2, counter, undom_verts,
                                            num_undom, grb_env, num_added,
                                            dom_sets, counter_cost);
                    // remove vertex from broadcast
                    rem_broadcast(vrow, vcol, 2, grid, dominated, v_ptr->len[1], v_ptr->dist[1]);   
                    // update the num_added and dom_sets
                    rem_from_dom_set_and_num_added(vcol, num_added);
                }
            }
        }
    }
    // try adding a broadcasting at strength 1
    // note we consider [0] since we are adding strength 1
    for (v_id = 0; v_id < ptr->num_dom_by[0]; v_id++) {
        // define coordinates of the vertex we want to add 
        // which dominates (row, col) 
        vrow = ptr->dom_by[0][v_id][0], vcol = ptr->dom_by[0][v_id][1];
        // if vertex is not already broadcasting
        if (grid[vrow][vcol] == 0) {
            v_ptr = &graph[0][vrow][vcol];
            // if vertex does not create a forbidden broadcast
            if (NoLocalForbiddenBroadcastStrength1(grid, *v_ptr)) { 
                // if there is not overlap in the broadcast by this vertex
                if (NoDominationOverlap(vcol, *v_ptr, 1, num_added, dom_sets)) {
                    // add this vertex to the broadcast
                    add_broadcast(vrow, vcol, 1, grid, dominated, v_ptr->len[0], v_ptr->dist[0]);
                    // update the num_added and dom_sets
                    add_to_dom_set_and_num_added(vcol, *v_ptr, 1, num_added, dom_sets);
                    // next step of recursion
                    ProduceAllSubBroadcasts(level + 1, grid, dominated, graph,
                                            cost + 1, counter, undom_verts,
                                            num_undom, grb_env, num_added,
                                            dom_sets, counter_cost);
                    // update the num_added and dom_sets
                    rem_broadcast(vrow, vcol, 1, grid, dominated, v_ptr->len[0], v_ptr->dist[0]);
                    // update the num_added and dom_sets
                    rem_from_dom_set_and_num_added(vcol, num_added);
                }
            }
        }
    }
    // vertex must be dominated
    return 0;
};

// test whether or not we find a contradiction for every subcase
bool ContradictionForEverySubCase(int grid[NUM_R][NUM_C], int dominated[NUM_R][NUM_C],
                                    struct Vertex graph[NUM_C][NUM_R][NUM_C], int nec_cost, 
                                    struct Counter counter[MAX_C + 1], GRBEnv grb_env,  
                                    int counter_cost) {
    int row, undom_verts[2*NUM_R][2];
    int num_undom = 0;

    // num_added[0] = number of vertices which are added to the broadcast to dominate column 4
    // num_added[1] = number of vertices which are added to the broadcast to dominate column NUM_C-5
    int num_added[2] = {0};
    // dom_sets[0][i] = bit-packed integers corresonding to the rows of the vertices dominated
    // in column 4 by the ith vertex added to dominate column 4
    // dom_sets[1][i] = bit-packed integers corresonding to the rows of the vertices dominated
    // in column NUM_C-5 by the ith vertex added to dominate column NUM_C-5
    uint32_t dom_sets[2][NUM_R];

    // determine the vertices undominated in columns 4 and NUM_C-5
    for (row = 0; row < NUM_R; row++) {
        if (dominated[row][4] == 0) {
            undom_verts[num_undom][0] = row;
            undom_verts[num_undom][1] = 4;
            num_undom += 1;
        }
        if (dominated[row][NUM_C-5] == 0) {
            undom_verts[num_undom][0] = row;
            undom_verts[num_undom][1] = NUM_C-5;
            num_undom += 1;
        }
    }
    // consider all possible subbroadcasts which extend this broadcast 
    // to dominate these remaining undominated vertices
    ProduceAllSubBroadcasts(0, grid, dominated, graph, nec_cost,
                            counter, undom_verts, num_undom, grb_env,
                            num_added, dom_sets, counter_cost);
    // if all contradiction were found
    if (counter[0].Fail == 0) {
        return 1;
    }
    return 0;
};
//
//****** CONSIDER ALL POSSIBLE SUBBROADCASTS ******//


//****** TRY TO FIND A CONTRADICTION ******//
//
bool ManageLPs(int grid[NUM_R][NUM_C], int dominated[NUM_R][NUM_C],
                struct Vertex graph[NUM_C][NUM_R][NUM_C], int cost,
                struct Counter counter[MAX_C + 1], GRBEnv grb_env) {
    int row;    
    int nec_cost;

    // first we test whether or not the vertices dominated by the base
    // broadcast can be dominated with less cost
    // original broadcast has cost = cost
    if (HasBroadcast(0, NULL, dominated, graph[0], cost - 1, grb_env)) {
        //****** FIGURES ******//
        if (PROD_FIG) {
            printf("col_del = []\n");
            printf("num_c = %d\n", NUM_C);
            printf("title = 'HasBroadcast'\n");
            ProduceGraphic();
        }
        //****** FIGURES ******//
        //
        // we found a conradiction
        counter[cost].HasBroadcast += 1;
        return 1;
    }

    // the vertices dominated by the base broadcast does not yield a
    // contradiction as is, try induciton
    if (InductiveArgument(dominated, graph, cost, grb_env)) {
        //****** FIGURES ******//
        if (PROD_FIG) {
            printf("title = 'InductiveArgument'\n");
            ProduceDominated(dominated);
            ProduceGraphic();
        }
        //****** FIGURES ******//
        //
        // we found a conradiction
        counter[cost].InductiveArgument += 1;
        return 1;
    }

    // unable to find a contradiction, add neccessary broadcast vertices
    nec_cost = AddNecBroadcast(grid, dominated, graph[0], cost);
    //****** FIGURES ******//
    if (PROD_FIG) {
        ProduceGrid(grid, nec_cost);
        printf("dominated = []\n");
        printf("title = 'NecessaryBroadcasts'\n");
        printf("col_del = []\n");
        printf("num_c = %d\n", NUM_C);
        ProduceGraphic();
        ProduceDominated(dominated);
    }
    //****** FIGURES ******//

    // first we test whether or not the vertices dominated by the 
    // broadcast with neccessary broadcast vertices can be dominated with less cost
    if (HasBroadcast(0, NULL, dominated, graph[0], nec_cost - 1, grb_env)) {
        //****** FIGURES ******//
        if (PROD_FIG) {
            printf("col_del = []\n");
            printf("num_c = %d\n", NUM_C);
            printf("title = 'NecessaryBroadcastsContradiction+HasBroadcast'\n");
            ProduceGraphic();
        }
        //****** FIGURES ******//
        //
        // remove neccessary broadcasts
        RemoveNecBroadcast(grid, dominated, graph[0]);
        // we found a conradiction
        counter[cost].NecessaryBroadcastHasBroadcast += 1;
        return 1;
    }

    // the vertices dominated by the broadcast with neccessary broadcast vertices
    // does not yield a contradiction as is, try induciton
    if (InductiveArgument(dominated, graph, nec_cost, grb_env)) {
        //****** FIGURES ******//
        if (PROD_FIG) {
            printf("title = 'NecessaryBroadcastsContradiction+InductiveArgument'\n");
            ProduceDominated(dominated);
            ProduceGraphic();
        }
        //****** FIGURES ******//
        //
        // we found a conradiction
        RemoveNecBroadcast(grid, dominated, graph[0]);
        counter[cost].NecessaryBroadcastInductiveArgument += 1;
        return 1;
    }
    
    // failed to prove anything with the neccessary broadcast vertices
    // now try to disprove all possible subcases
    if (ContradictionForEverySubCase(grid, dominated, graph, nec_cost,
                                     counter, grb_env, cost)) {
        RemoveNecBroadcast(grid, dominated, graph[0]);
        // we found a contradiction
        return 1;
    }

    // unable to find a contradiction
    // remove neccessary broadcasts
    RemoveNecBroadcast(grid, dominated, graph[0]);
    return 0;
};
//
//****** TRY TO FIND A CONTRADICTION ******//


//****** TOOLS FOR CREATING ALL THE BROADCASTS ******//
//
// test whether given broadcast is cannonically greatest given the 
// group actions on the Cartesian product of two paths
bool PxP_Not_Cannonical(int grid[NUM_R][NUM_C]) {
    int row, col;
    // horizontal flip
    for (row = 0; row < NUM_R; row++) {
        for (col = 0; col < NUM_C; col++) {
            // current grid is cannonically greater
            if (grid[row][col] > grid[row][NUM_C - 1 - col]) {
                // break the loops
                goto End1;
            }
            // not cannonically greater
            else if (grid[row][col] < grid[row][NUM_C - 1 - col]) {
                return 1;
            }
        }
    }
    End1:
    // vertical flip
    for (row = 0; row < NUM_R; row++) {
        // only dealing with non-padded columns
        for (col = 4; col < NUM_C - 4; col++) {
            // current grid is cannonically greater
            if (grid[row][col] > grid[NUM_R - 1 - row][col]) {
                // break the loops
                goto End2;
            }
            // not cannonically greater
            else if (grid[row][col] < grid[NUM_R - 1 - row][col]) {
                return 1;
            }
        }
    }
    End2:
    // horizontal and vertical flip
    for (row = 0; row < NUM_R; row++) {
        // only dealing with non-padded columns
        for (col = 4; col < NUM_C - 4; col++) {
            // current grid is cannonically greater
            if (grid[row][col] > grid[NUM_R - 1 - row][NUM_C - 1 - col]) {
                // break the loops
                goto End3;
            }
            // not cannonically greater
            else if (grid[row][col] < grid[NUM_R - 1 - row][NUM_C - 1 - col]) {
                return 1;
            }
        }
    }
    End3:
    // made it through the checks, grid is cannonically greatest form
    return 0;
};

// check whether a given broadcast dominates the middle of the grid
bool DoesNotDominate(int dominated[NUM_R][NUM_C]) {
    int row, col;

    // center columns need to be dominated
    for (row = 0; row < NUM_R; row++) {
        for (col = 6; col < NUM_C - 6; col++) {
            // not dominated 
            if (dominated[row][col] == 0) {
                return 1;
            }
        }
    }
    return 0;
};

// check if there is a forbidden broadcast
bool ForbiddenBroadcast(int grid[NUM_R][NUM_C],
                        struct Vertex graph[NUM_R][NUM_C]) {
    int row, col, i;
    struct Vertex *graph_ptr;

    for (row = 0; row < NUM_R; row++) {
        for (col = 0; col < NUM_C; col++) {
            // if vertex is broadcasting at strength 1. Note that we need not
            // check for those of strength 2 since they will be caught by the 
            // strength 1 ones
            if (grid[row][col] == 1) {
                graph_ptr = &graph[row][col];
                // if it has a neighbor broadcasting at non-zero strength
                for (i = 0; i < graph_ptr->len_E[0]; i++) {
                    if (grid[graph_ptr->dist_E[0][i][0]][graph_ptr->dist_E[0][i][1]]!= 0) {
                        return 1;
                    }
                }
                // if it has a vertex at distance 2 broadcasting at strength 1
                for (i = 0; i < graph_ptr->len_E[1]; i++) {
                    if (grid[graph_ptr->dist_E[1][i][0]][graph_ptr->dist_E[1][i][1]] == 1) {
                        return 1;
                    }
                }
            }
        }
    }
    // no forbidden broadcast found
    return 0;
};
//****** TOOLS FOR CREATING ALL THE BROADCASTS ******//


//****** RECURSION TO CONSTRUCT BROADCAST ******//
//
bool fill_grid(int level, int grid[NUM_R][NUM_C],
              int dominated[NUM_R][NUM_C], int cost, int des_cost,
              struct Vertex graph[NUM_C][NUM_R][NUM_C], 
              struct Counter counter[MAX_C + 1], GRBEnv grb_env) {
    int row, col;
    struct Vertex *graph_ptr;

    // if we failed to prove a case
    if (counter[0].Fail != 0) {
        return 0;
    }

    // define row and column of current vertex
    // row is floor of division of level by
    // number of columns less the 8 padded columns
    row = level/(NUM_C - 8);
    // col corresonds with level mod (number of columns less
    // the 8 padded columns) + 4 padded columns to shift back over
    col = (level % (NUM_C - 8)) + 4;

    // if broadcast is of desired cost
    if (des_cost == cost) {
        // is grid not cannonically greatest?
        if (PxP_Not_Cannonical(grid)) {
            return 0;
        }
        // otherwise cannonical
        counter[des_cost].C += 1;
        // does the broadcast dominate the middle of the grid?
        if (DoesNotDominate(dominated)) {
            counter[des_cost].DoesNotDominate += 1;
            return 0;
        }
        // does it contain a forbidden broadcast?
        if (ForbiddenBroadcast(grid, graph[0])) {
            counter[des_cost].ForbiddenBroadcast += 1;
            return 0;
        }
        //****** FIGURES ******//
        if (PROD_FIG) {
            ProduceGrid(grid, cost);
            printf("dominated = []\n");
            printf("title = 'BaseCase'\n");
            printf("col_del = []\n");
            printf("num_c = %d\n", NUM_C);
            ProduceGraphic();
            ProduceDominated(dominated);
        }
        //****** FIGURES ******//
        
        // hand to the contradiction stage
        if (ManageLPs(grid, dominated, graph, cost, counter,
                        grb_env)) {
            return 0;
        } else {
            counter[0].Fail += 1;
            ProduceGrid(grid, cost);
            printf("dominated = []\n");
            printf("title = 'Fail'\n");
            printf("col_del = []\n");
            printf("num_c = %d\n", NUM_C);
            ProduceGraphic();
            return 0;
        }
    }

    // if ran out of rows
    if (row == NUM_R) {
        return 0;
    }

    graph_ptr = &graph[0][row][col];

    // TRY BROADCASTING AT STRENGTH 2
    if (des_cost - cost >= 2) {
        // add strength 2 vertex
        add_broadcast(row, col, 2, grid, dominated, graph_ptr->len[1],
                        graph_ptr->dist[1]);
        // move to next vertex
        fill_grid(level + 1, grid, dominated, cost + 2,
                    des_cost, graph, counter, grb_env);
        // remove strength 2 vertex
        rem_broadcast(row, col, 2, grid, dominated, graph_ptr->len[1],
                        graph_ptr->dist[1]);    
    }

    // TRY BROADCASTING AT STRENGTH 1
    // add strength 1 vertex
    add_broadcast(row, col, 1, grid, dominated, graph_ptr->len[0],
                    graph_ptr->dist[0]);
    // move to next vertex
    fill_grid(level + 1, grid, dominated, cost + 1,
                des_cost, graph, counter, grb_env);
    // remove strength 2 vertex
    rem_broadcast(row, col, 1, grid, dominated, graph_ptr->len[0],
                    graph_ptr->dist[0]);    

    // TRY BROADCASTING AT STRENGTH 0
    fill_grid(level + 1, grid, dominated, cost,
                des_cost, graph, counter, grb_env);

    // no other option
    return 0;
};
//
//****** RECURSION TO CONSTRUCT BROADCAST ******//



//****** SET EVERYTHING RUNNING ******//
//
int main(void) {
    int cost;
    // array containing broadcast, [x][y] = i iff vertex at row x and col y is broadcasting
    // at strength i
    int grid[NUM_R][NUM_C] = { 0 };
    // array to indicate which vertices are dominated, [x][y] != 0 iif vertex at row x
    // and col y is dominated
    int dominated[NUM_R][NUM_C] = { 0 };
    // array containing the graph information, indexed by the number of columns deleted
    // from the graph, ROWS, COLS, then a struct containing the relevant information
    // for each vertex
    struct Vertex graph[NUM_C][NUM_R][NUM_C] = {};

    //****** FIGURES ******//
    if (PROD_FIG) {
        printf("from matplotlib.backends.backend_pdf import PdfPages\n");
        printf("from graphic import produce_graphic\n");
        printf("#");
    }
    //****** FIGURES ******//
    
    // define a gurobi enviroment 
    GRBEnv grb_env = GRBEnv();
    grb_env.set(GRB_IntParam_OutputFlag,0);

    // fill in graph information
    make_graph(graph);

    // make counter
    struct Counter counter[MAX_C + 1] = {};

    //****** FIGURES ******//
    if (PROD_FIG) {
        if (CYCLE == 1) {
            printf("pdf = PdfPages('C%d x C%d, cost = %d -- %d.pdf')\n", NUM_R, NUM_C, MIN_C, MAX_C);
        } else if (CYCLE == 0) {
            printf("pdf = PdfPages('P%d x C%d, cost = %d -- %d.pdf')\n", NUM_R, NUM_C, MIN_C, MAX_C);
        } else {
            printf("Check CYCLE val");
            return 0;
        }
    } else {
        if (CYCLE == 1) {
            printf("C%d x C%d, cost = %d -- %d\n", NUM_R, NUM_C, MIN_C, MAX_C);
        } else if (CYCLE == 0) {
            printf("'P%d x C%d, cost = %d -- %d\n", NUM_R, NUM_C, MIN_C, MAX_C);
        } else {
            printf("Check CYCLE val");
            return 0;
        }
    }
    //****** FIGURES ******//

    // call recursion
    for (cost = MIN_C; cost <= MAX_C; cost ++) {
        fill_grid(0, grid, dominated, 0, cost, graph, counter, grb_env);
        // output results
        if (PROD_FIG) {
            printf("#### Cost = %d ###\n", cost);
            printf("#|C|: %llu\n", counter[cost].C);
            printf("#DoesNotDominate: %llu\n", counter[cost].DoesNotDominate);
            printf("#ForbiddenBroadcast: %llu\n", counter[cost].ForbiddenBroadcast);
            printf("#HasBroadcast: %llu\n", counter[cost].HasBroadcast);
            printf("#InductiveArgument: %llu\n", counter[cost].InductiveArgument);
            printf("#NecessaryBroadcastHasBroadcast: %llu\n", counter[cost].NecessaryBroadcastHasBroadcast);
            printf("#NecessaryBroadcastInductiveArgument: %llu\n", counter[cost].NecessaryBroadcastInductiveArgument);
            printf("#AllSubcasesHasBroadcast: %llu\n", counter[cost].AllSubcasesHasBroadcast);
            printf("#AllSubcasesInductiveArgument: %llu\n", counter[cost].AllSubcasesInductiveArgument);
            printf("#Failed: %d\n",counter[cost].Fail);
        } else {
            printf("### Cost = %d ###\n", cost);
            printf("|C|: %llu\n", counter[cost].C);
            printf("DoesNotDominate: %llu\n", counter[cost].DoesNotDominate);
            printf("ForbiddenBroadcast: %llu\n", counter[cost].ForbiddenBroadcast);
            printf("HasBroadcast: %llu\n", counter[cost].HasBroadcast);
            printf("InductiveArgument: %llu\n", counter[cost].InductiveArgument);
            printf("NecessaryBroadcastHasBroadcast: %llu\n", counter[cost].NecessaryBroadcastHasBroadcast);
            printf("NecessaryBroadcastInductiveArgument: %llu\n", counter[cost].NecessaryBroadcastInductiveArgument);
            printf("AllSubcasesHasBroadcast: %llu\n", counter[cost].AllSubcasesHasBroadcast);
            printf("AllSubcasesInductiveArgument: %llu\n", counter[cost].AllSubcasesInductiveArgument);
            printf("Failed: %d\n",counter[cost].Fail);
        }
    }

    if (PROD_FIG) {
        if (counter[0].Fail == 0) {
            printf("#### Proven ####\n");
        } else {
            printf("#### Failed to prove ####\n");
        }
    } else {
        if (counter[0].Fail == 0) {
            printf("#### Proven ####\n");
        } else {
            printf("#### Failed to prove ####\n");
        }
    }

    //
    //****** FIGURES ******//
    if (PROD_FIG) {
        printf("pdf.close()");
    }
    //****** FIGURES ******//
    
};
//
//****** SET EVERYTHING RUNNING ******//
