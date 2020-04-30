#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <stack>

using namespace std;

class Dimension {
private:
    double *dim;

public:
    explicit Dimension(int dimensions) {
        dim = new double[dimensions];
    }

    double getNumber(int pos) {
        return dim[pos];
    }

    void setNumber(int pos, double num) {
        dim[pos] = num;
    }
};

class Center {
private:
    double *v;

public:
    explicit Center(int centers) {
        v = new double[centers];
    }

    double getNumber(int pos) {
        return v[pos];
    }

    void setNumber(int pos, double num) {
        v[pos] = num;
    }
};

//Boundary data structure for dimension axis
struct BND {
    double max;
    double min;
};

//Get data from file location and put it into flexible vector
vector<Dimension *> getData(int start_point, int count, int dim, char file_loc[]) {
    vector<Dimension *> d;

    d.reserve(dim);
    for (int i = 0; i < dim; i++) {
        d.push_back(new Dimension(count));
    }
    ifstream fin(file_loc);

    //Universal way for extracting data(Only dependence is knowing dimension inside file)
    double num;
    int it = 0;//iterations
    int c_dim = 0; //current dimension
    for (int i = 0; fin >> num; i++) {
        if (i < start_point * (dim + 1) - dim) continue; //skip until start point
        if ((i == 0 || i % (dim + 1) == 0) && i != 1)
            continue; // skip if it is initials, and not the first row
        if (it == count) break;
        d.at(c_dim)->setNumber(it, num);
        if (c_dim == dim - 1) {
            c_dim = 0;
            it++;
            continue;
        }
        c_dim++;
    }
    return d;
}

vector<BND> getBndr(int start_point, int count, int dim, char file_loc[]) {
    vector<BND> boundaries;
    boundaries.reserve(dim);
    for (int k = 0; k < dim; k++) {
        boundaries.push_back(BND{0, 0});
    }
    vector<Dimension *> dimension = getData(start_point, count, dim, file_loc);
    double cmax, cmin;
    for (int i = 0; i < dim; i++) {
        cmax = dimension.at(i)->getNumber(0);
        cmin = dimension.at(i)->getNumber(0);
        for (int j = 0; j < count; j++) {
            if (cmax < dimension.at(i)->getNumber(j)) {
                cmax = dimension.at(i)->getNumber(j);
            }
            if (cmin > dimension.at(i)->getNumber(j)) {
                cmin = dimension.at(i)->getNumber(j);
            }
        }
        boundaries.at(i).max = cmax;
        boundaries.at(i).min = cmin;
    }
    return boundaries;
}

stack<double>
getDensity(char *file, int num_points, int start_point, vector<BND> bndr, int num_divider, int dim) {
    //points
    vector<Dimension *> dimension = getData(start_point, num_points, dim, file);

    //cluster centers
    double step;//we use equally spaced cluster centers - grid cells for each axis
    double sp, eps = 0.001;
    int num_cells = pow(num_divider + 1, dim);

    int M = 2; //fuzzifier
    vector<Center *> v;
    v.reserve(dim);
    for (int k = 0; k < dim; k++) {
        v.push_back(new Center(num_cells));
    }
    int changer = 0;

    //Grid cell maker based on dimension
    //trigger change every pow(num_divider + 1, i), but max change num_divider + 1
    for (int i = 0; i < dim; i++) {
        sp = bndr.at(i).min;
        step = (bndr.at(i).max - bndr.at(i).min) / num_divider;
        for (int j = 0; j < num_cells; j++) {
            changer = pow(num_divider + 1, i);

            if (j % changer == 0 && sp >= bndr.at(i).max - eps)
                sp = bndr.at(i).min;
            else if (j % changer == 0 && j != 0)
                sp += step;
            v.at(i)->setNumber(j, sp);
        }
    }

    //Membership matrix - Fuzzy Partition Matrix
    double FM[num_points][num_cells];
    double nom, den = 0.0, tDen;
    for (int i = 0; i < num_points; i++) {
        for (int j = 0; j < num_cells; j++) {
            tDen = 0.0;

            //sum of the reciprocal of denominator of denominator
            //for clusters
            for (int m = 0; m < num_cells; m++) {
                //sum of the squared parts in denominator of denominator
                for (int k = 0; k < dim; k++) {
                    den += pow(abs((dimension.at(k)->getNumber(i) - v.at(k)->getNumber(m))), 2 / (M - 1));
                }
                tDen += 1 / den;
                den = 0.0;
            }
            nom = 0.0;
            //sum of the squared parts in nominator of denominator
            for (int k = 0; k < dim; k++) {
                nom += pow(abs((dimension.at(k)->getNumber(i) - v.at(k)->getNumber(j))), 2 / (M - 1));
            }

            FM[i][j] = 1 / (nom * tDen);
            if (isnan(FM[i][j])) FM[i][j] = 1; // check if it has infinite membership degree or, nom, tDen is 0
        }
    }

    //Initialise and fill density stack
    stack<double> density;
    double sum = 0;
    for (int i = 0; i < num_cells; i++) {
        for (int j = 0; j < num_points; j++) {
            sum += FM[j][i];
        }
        density.push(sum);
        sum = 0;
    }
    return density;
}

int main() {
    char _file[50];
    int _num_points, _dim, _start_point_DCM, _start_point_WCM, _div;

    cout << "Please enter the file location: ";
    cin >> _file;

    cout << "Please enter number of points: ";
    cin >> _num_points;

    cout << "Please enter number of dimension: ";
    cin >> _dim;

    _start_point_DCM = 630; //from which element you want to "cut" window
    _start_point_WCM = 1; //from which element you want to "cut" window

    //get boundaries for each axis(dimensions)
    vector<BND> bndr = getBndr(_start_point_DCM, _num_points, _dim, _file);

    //Below declarations are for creating grid cells
    _div = 4; // divider of each axis(More divider, more accuracy but slower speed)

    clock_t tStart = clock();
    stack<double> dcm = getDensity(_file, _num_points, _start_point_DCM, bndr, _div, _dim);
    stack<double> wcm;

    double sum = 0;
    int incr_step = 1; //increment size for windows
    for (int w = 0; w < 1000; w += incr_step) { //windows count
        wcm = getDensity(_file, _num_points, _start_point_WCM + w, bndr, _div, _dim); //update window

        stack<double> s_temp = dcm;
        for (int i = 0; i < pow(_div + 1, _dim); i++) {
            sum += abs(s_temp.top() - wcm.top());
            s_temp.pop();
            wcm.pop();
        }

        double SL = 1.0 - sum / (2 * _num_points);
        cout << w + 1 << ". window: " << SL << endl;

        sum = 0;
    }
    printf("Time taken: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);
    return 0;
}
