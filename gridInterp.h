#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

namespace GRID_INTERP{
    /**
     * @brief There are two types of interpolation modes.
     * In POINT MODE the interpolated values coinside with the coordinates. 
     * For every coordinate there has to be one value. 
     * Therefore the number of values is equal to the number of coordinates
     * 
     * In CELL MODE the interpolated values correspond to cells and the coordinates to 
     * the interface between the cells. In CELL MODE the coordinates are one plus the number of values.
     * eg.  | v1 | v2 | ... | vn-1 | vn |
     *      X1   X2  X3    xn-1   xn    xn+1
     * 
     */
    enum MODE {CELL, POINT, INVALID_MODE};

    /**
     * @brief There are two interpolation methods.
     * Linear and nearest neighborhood.
     * 
     */
    enum METHOD{LINEAR, NEAREST, INVALID_METHOD};

    /**
     * @brief The interpolation type it is used only in 3D. For 1 and 2 D it is ignored
     * The normal type is when interpolation is using a full 3D grid.
     * The layer type is used when the elevation of the 3D field is variable and cannot be defined as one value
     * 
     */
    enum TYPE{GRID, LAYER, INVALID_TYPE};

    /**
     * @brief An axis class is a container for axis tick values.
     * 
     * The tick values can have variable space or constant.
     * If the space between the tick values is constant it is advantageous 
     * to use the setAxis for constant ticks. 
     * 
     * An axis object is used to describe the location of the values during interpolation
     * therefore the number of ticks must be consistent with the grid values.
     * For the GRID_INTERP::MODE::POINT mode the number of ticks must be exaclty the number of values 
     * For the GRID_INTERP::MODE::CELL mode the number of ticks is the number of values + 1
     * 
     */
    class axis{
    public:
        //! The constructor does absolutely nothing
        axis(){};
        /**
         * @brief Set the Axis tick values when the space between them is constant 
         * 
         * @param origin is the starting value 
         * @param dx is the space between the ticks
         * @param N is the number of tick values to generate.
         */
        void setAxis(double origin, double dx, int N);
        /**
         * @brief Set the Axis tick values with the input vector.
         * This is used to set variable space grid data
         * 
         * @param x_in 
         */
        void setAxis(std::vector<double> &x_in);
        /**
         * @brief Finds the lower and upper ids that enclose the input value x_in.
         * Also returns the parametric value that correspons to x1
         * 
         * @param x_in is the query point
         * @param i1 is the id of the lower limit
         * @param i2 is the id of the upper limit
         * @param t is the parametric value so that x_in = x[i1]*(1-t) + x[i2]*t
         */
        void findIIT(double x_in, int &i1, int &i2, double &t);
        //! so far this is not used so I might delete this
        void findIcell(double x_in, int &i);
        /**
         * @brief Reads the axis data from file
         * 
         * The data are organized in the following format:
         * 
         * 1st line: DATA_FORMAT N
         * 
         * where DATA_FORMAT is either CONST or VAR and N is the number of ticks
         * 
         * If the DATA_FORMAT is CONST the second line consists of 2 numbers
         * origin dx
         * If the DATA_ORMAT is VAR then we read N lines where each line contains the values of the ticks
         * 
         * @param filename 
         */
        void dataFromFile(std::string filename);
        //! Delets all data and resest the axis object to its original state
        void reset();
        //! accessor for the i element of the axis
        double& operator()(int i);
        //! returns the index that access the last element of the tick values
        //! If ax is an axis object then ax(ax.last()) is equal to ax.x[ax.x.size()-1]
        int last();

    private:
        //! is the vector that containts the tick values
        std::vector<double> x;
        //! is the index of the last element
        int N = 0;
        //! this is the origin of the grid. When the object is populated with variable data this is not used
        double origin = 0.0;
        //! The distance between the ticks. This is used only for constant grid values
        double dx = 0.0;
        //! This flag is set to true or false according to which setAxis method is used.
        bool bConst = false;
        //! Finds cell searches the segment that containts the x. Is doing so by dividing into half and repeating it self
        void findCell(int &i, int &ii, double x);
    };

    void axis::reset(){
        x.clear();
        origin = 0.0;
        dx = 0.0;
        bConst = false;
        N = 0;
    }
    double& axis::operator()(int i){
        return x[i];
    }
    int axis::last(){
        return x.size()-1;
    }

    void axis::setAxis(double origin_in, double dx_in, int n){
        x.clear();
        for(int i = 0; i < n; i++){
            x.push_back(origin_in + static_cast<double>(i)*dx_in);
        }
        bConst = true;
        origin = origin_in;
        dx = dx_in;
        N = x.size() - 1;
    }

    void axis::setAxis(std::vector<double> &x_in){
        x = x_in;
        N = x.size() - 1;
    }

    void axis::findCell(int &i, int &ii, double x_in){
        switch (bConst)
        {
        case true:
            {
                i = static_cast<int>((x_in - origin)/dx);
                ii = i + 1;
            }
            break;
        default:
            {
                int mid = (i + ii)/2;
                if (x_in <= x[mid])
                    ii = mid;
                else
                    i = mid;
                
                if (ii - i <= 1)
                    return;
                else
                    findCell(i, ii, x_in);
            }
            break;
        }
    }

    void axis::findIIT(double x_in, int &i1, int &i2, double &t){
        if (x_in < x[0]){
            i1 = 0;
            i2 = 0;
            t = 0.0;
            return;
        }
        if (x_in > x[N]){
            i1 = N;
            i2 = N;
            t = 1.0;
            return;
        }
        i1  = 0;
        i2 = N;
        findCell(i1, i2, x_in);
        t = (x_in - x[i1])/(x[i2] - x[i1]);
    }

    void axis::findIcell(double x_in, int &i){
        if (x_in < x[0]){
            i = 0;
            return;
        }
        if (x_in > x[x.size()-1]){
            i = N;
            return;
        }
        int i1  = 0;
        int i2 = N;
        findCell(i1, i2, x_in);
        i = i1;
    }

    void axis::dataFromFile(std::string filename){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file" << filename << std::endl;
        }
        else{
            int Nticks;
            std::string line, tmp;
            getline(datafile, line);
            {
                std::istringstream inp(line.c_str());
                inp >> tmp;
                
                if (tmp.compare("CONST") == 0)
                    bConst = true;
                else
                    bConst = false;

                inp >> Nticks;
            }

            if (bConst){
                getline(datafile, line);
                std::istringstream inp(line.c_str());
                inp >> origin;
                inp >> dx;
                setAxis(origin, dx, Nticks);
            }
            else{
                x.clear();
                double xx;
                for (int i = 0; i < Nticks; i++){
                    getline(datafile, line);
                    std::istringstream inp(line.c_str());
                    inp >> xx;
                    x.push_back(xx);
                }
                N = x.size() - 1;
            }
        }
    }



    class GridValues{
    public:
        GridValues(){};
        double& operator()(int l, int r, int c);
        double& operator()(int r, int c);
        double& operator()(int c);
        void resize(int nl, int nr, int nc);
        void set(int l, int r, int c, double v_in);
        void reset();
        int nx();
        int ny();
        int nz();
    private:
        std::vector<std::vector<std::vector<double>>> v;
    };

    double& GridValues::operator()(int l, int r, int c){
        return v[l][r][c];
    }
    double& GridValues::operator()(int r, int c){
        return v[0][r][c];
    }
    double& GridValues::operator()(int c){
        return v[0][0][c];
    }
    void GridValues::resize(int nl, int nr, int nc){
        v.clear();
        //std::vector<double>(nc,0);
        //std::vector<std::vector<double>>(nr, std::vector<double>(nc,0));
        v.resize(nl,std::vector<std::vector<double>>(nr, std::vector<double>(nc,0)));
    }
    void GridValues::set(int l, int r, int c, double v_in){
        v[l][r][c] = v_in;
    }

    void GridValues::reset(){
        v.clear();
    }

    int GridValues::nx(){
        return v[0][0].size();
    }
    int GridValues::ny(){
        return v[0].size();
    }
    int GridValues::nz(){
        return v.size();
    }

    template<int dim>
    class interp{
    public:
        interp();
        double interpolate(double x, double y = 0, double z = 0);
        void setAxis(int idim, std::vector<double> &x_in);
        void setAxis(int idim, double origin_in, double dx_in, int n);
        void setModeMethod(MODE mode_in, METHOD method_in);
        void getDataFromFile(std::string filename);
        void reset();

    private:
        std::vector<axis> a;
        GridValues v;
        GridValues elev;
        MODE mode = MODE::INVALID_MODE;
        METHOD method = METHOD::INVALID_METHOD;
        double interp1D(double x);
        double interp2D(double x, double y);
        double interp3D(double x, double y, double z);

        void cellLinearcorrect(int idim, double x, int &i1, int &i2, double &x1, double &x2);
    };

    template<int dim>
    interp<dim>::interp(){
        a.clear();
        a.resize(dim);
    }

    template<int dim>
    void interp<dim>::reset(){
        for(int idim = 0; idim < dim; ++ idim)
            a[idim].reset();
        v.reset();
        mode = MODE::INVALID_MODE;
        method = METHOD::INVALID_METHOD;
    }

    template <int dim>
    void interp<dim>::setAxis(int idim, std::vector<double> &x_in){
        if (idim < dim)
            a[idim].setAxis(x_in);
    }

    template <int dim>
    void interp<dim>::setAxis(int idim, double origin_in, double dx_in, int n){
        if (idim < dim)
            a[idim].setAxis(origin_in, dx_in, n);
    }

    template<int dim>
    void interp<dim>::setModeMethod(MODE mode_in, METHOD method_in){
        mode = mode_in;
        method = method_in;
    }

    template <int dim>
    void interp<dim>::getDataFromFile(std::string filename){
        std::ifstream datafile(filename.c_str());
        if (!datafile.good()) {
            std::cout << "Can't open the file" << filename << std::endl;
        }
        else{
            MODE md;
            METHOD mthd;
            std::string line, tmp;
            std::vector<int> dims;
            
            getline(datafile, line);
            {// MODE METHOD NX NY NZ
                std::istringstream inp(line.c_str());
                // Get mode
                inp >> tmp;
                if (tmp.compare("POINT") == 0)
                    md = MODE::POINT;
                else if (tmp.compare("CELL") == 0)
                    md = MODE::CELL;
                else
                    md = MODE::INVALID_MODE;

                //get method
                inp >> tmp;
                if (tmp.compare("LINEAR") == 0)
                    mthd = METHOD::LINEAR;
                else if (tmp.compare("NEAREST") == 0)
                    mthd = METHOD::NEAREST;
                else
                    mthd = METHOD::INVALID_METHOD; 
                int d;
                for(int idim = 0; idim < dim; ++idim){
                    inp >> d;
                    dims.push_back(d);
                }
                if (dim == 1)
                    v.resize(1,1,dims[0]);
                else if (dim == 2)
                    v.resize(1,dims[1],dims[0]);
                else if (dim == 3)
                    v.resize(dims[2], dims[1], dims[0]);
            }

            getline(datafile, line);
            {// FILEX, FILEY, FILEZ
                std::istringstream inp(line.c_str());
                for(int i = 0; i < dim; ++i){
                    inp >> tmp;
                    a[i].dataFromFile(tmp);
                }
            }
            if (dim == 1){
                double vv;
                for (int i = 0; i < dims[0]; ++i){
                    getline(datafile, line);
                    std::istringstream inp(line.c_str());
                    inp >> vv;
                    v.set(0,0,i,vv);
                }
            }
            else{
                int nz = 1;
                int ny, nx;
                if (dim == 2){
                    nx = dims[0];
                    ny = dims[1];
                }
                else if (dim == 3){
                    nx = dims[0];
                    ny = dims[1];
                    nz = dims[2];

                }
                double vv;
                for (int k = 0; k < nz; ++k)
                    for (int i = 0; i < ny; ++i){
                        getline(datafile, line);
                        std::istringstream inp(line.c_str());
                        for (int j = 0; j < nx; ++j){
                            inp >> vv;
                            v.set(k, i, j, vv);
                        }
                    }
            }
            setModeMethod(md, mthd);
        }
    }

    template<int dim>
    double interp<dim>::interpolate(double x, double y, double z){
        double outcome;
        switch (dim)
        {
        case 1:
            outcome = interp1D(x);
            break;
        case 2:
            outcome = interp2D(x,y);
            break;
        case 3:
            outcome = interp3D(x,y,z);
            break;
        default:
            outcome = std::numeric_limits<double>::quiet_NaN();
            break;
        }
        return outcome;
    }

    template<int dim>
    double interp<dim>::interp1D(double x){
        double outcome;
        int i1, i2;
        double t;
        a[0].findIIT(x,i1, i2,t);
        
        switch (mode)
        {
        case MODE::POINT:

            switch (method)
            {
            case METHOD::LINEAR:
                outcome = v(i1)*(1-t) + v(i2)*t;
                break;
            case METHOD::NEAREST:
                if (std::abs(a[0](i1) - x) <= std::abs(a[0](i2) - x))
                    outcome = v(i1);
                else
                    outcome = v(i2);
                break;
            default:
                outcome = std::numeric_limits<double>::quiet_NaN();
                break;
            }
            break;
        case MODE::CELL:
            switch (method)
            {
            case METHOD::LINEAR:
                {
                    double x1, x2;
                    cellLinearcorrect(0, x, i1, i2, x1, x2);
                    double tt = (x - x1)/(x2 - x1);
                    outcome = v(i1)*(1-tt) + v(i2)*tt;
                }
                break;
            case METHOD::NEAREST:
                if (i1 >= v.nx())
                    i1 = v.nx() - 1;
                outcome = v(i1);
                break;
            default:
                outcome = std::numeric_limits<double>::quiet_NaN();
                break;
            }
            break;
        default:
            outcome = std::numeric_limits<double>::quiet_NaN();
            break;
        }
        return outcome;
    }
  

    template<int dim>
    double interp<dim>::interp2D(double x, double y){
        double outcome;
        int i1, i2, j1, j2;
        double tx, ty;
        a[0].findIIT(x, j1, j2, tx);
        a[1].findIIT(y, i1, i2, ty);
        switch (mode)
        {
        case MODE::POINT:
            {
                switch (method)
                {
                case METHOD::LINEAR:
                    {
                        double y1 = v(i1, j1)*(1-tx) + v(i1, j2)*tx;
                        double y2 = v(i2, j1)*(1-tx) + v(i2, j2)*tx;
                        outcome = y1*(1-ty) + y2*ty;
                        
                    }
                    break;
                case METHOD::NEAREST:
                    {
                        int i, j;
                        if (std::abs(a[0](j1) - x) <= std::abs(a[0](j2) - x))
                            j = j1;
                        else
                            j = j2;

                        if (std::abs(a[1](i1) - y) <= std::abs(a[1](i2) - y))
                            i = i1;
                        else
                            i = i2;
                        outcome = v(i,j);
                    }
                    break;
                default:
                    break;
                }
            }
            break;
        case MODE::CELL:
            {
                switch (method)
                {
                case METHOD::LINEAR:
                    {
                        double x1, x2, y1, y2;
                        cellLinearcorrect(0, x, j1, j2, x1, x2);
                        cellLinearcorrect(1, y, i1, i2, y1, y2);
                        
                        double tx = (x - x1)/(x2 - x1);
                        double ty = (y - y1)/(y2 - y1);

                        // y1 y2 correspond to values here
                        y1 = v(i1, j1)*(1-tx) + v(i1, j2)*tx;
                        y2 = v(i2, j1)*(1-tx) + v(i2, j2)*tx;
                        outcome = y1*(1-ty) + y2*ty;
                    }
                    break;
                case METHOD::NEAREST:
                    {
                        if (i1 >= v.ny())
                            i1 = v.ny() - 1;
                        if (j1 >= v.nx())
                            j1 = v.nx() - 1;
                        outcome = v(i1, j1);
                    }
                    break;
                default:
                    outcome = std::numeric_limits<double>::quiet_NaN();
                    break;
                }
            }
            break;
        default:
            outcome = std::numeric_limits<double>::quiet_NaN();
            break;
        }
        return outcome;
    }

    template<int dim>
    double interp<dim>::interp3D(double x, double y, double z){
        int i1, i2, j1, j2, k1, k2;
        double tx, ty, tz;
        a[0].findIIT(x, j1, j2, tx);
        a[1]. findIIT(y, i1, i2, ty);
        a[2]. findIIT(z, k1, k2, tz);
        switch (mode)
        {
        case MODE::POINT:
            {
                switch (method)
                {
                case METHOD::LINEAR:
                    {
                        double z1y1 = v(k1, i1, j1)*(1-tx) + v(k1, i1, j2)*tx;
                        double z1y2 = v(k1, i2, j1)*(1-tx) + v(k1, i2, j2)*tx;
                        double z1 = z1y1*(1-ty) + z1y2*ty;

                        double z2y1 = v(k2, i1, j1)*(1-tx) + v(k2, i1, j2)*tx;
                        double z2y2 = v(k2, i2, j1)*(1-tx) + v(k2, i2, j2)*tx;
                        double z2 = z2y1*(1-ty) + z2y2*ty;

                        return z1*(1-tz) + z2*tz;
                    }
                    break;
                case METHOD::NEAREST:
                    {
                        int i, j, k;
                        if (std::abs(a[0](j1) - x) <= std::abs(a[0](j2) - x))
                            j = j1;
                        else
                            j = j2;

                        if (std::abs(a[1](i1) - y) <= std::abs(a[1](i2) - y))
                            i = i1;
                        else
                            i = i2;
                        if (std::abs(a[2](k1) - y) <= std::abs(a[2](k2) - y))
                            k = k1;
                        else
                            k = k2;
                        return v(k, i, j);

                    }
                    break;
                default:
                    break;
                }
            }
            break;
        case MODE::CELL:
            {
                switch (method)
                {
                case METHOD::LINEAR:
                    {
                        double x1, x2, y1, y2, z1, z2;
                        cellLinearcorrect(0, x, j1, j2, x1, x2);
                        cellLinearcorrect(1, y, i1, i2, y1, y2);
                        cellLinearcorrect(2, z, k1, k2, z1, z2);
                        
                        double tx = (x - x1)/(x2 - x1);
                        double ty = (y - y1)/(y2 - y1);
                        double tz = (z - z1)/(z2 - z1);

                        // z1 and z2 correspond to values here
                        double z1y1 = v(k1, i1, j1)*(1-tx) + v(k1, i1, j2)*tx;
                        double z1y2 = v(k1, i2, j1)*(1-tx) + v(k1, i2, j2)*tx;
                        z1 = z1y1*(1-ty) + z1y2*ty;

                        double z2y1 = v(k2, i1, j1)*(1-tx) + v(k2, i1, j2)*tx;
                        double z2y2 = v(k2, i2, j1)*(1-tx) + v(k2, i2, j2)*tx;
                        z2 = z2y1*(1-ty) + z2y2*ty;

                        return z1*(1-tz) + z2*tz;
                    }
                    break;
                case METHOD::NEAREST:
                    {   
                        if (k1 >= v.nz())
                            k1 = v.nz() - 1;
                        if (i1 >= v.ny())
                            i1 = v.ny() - 1;
                        if (j1 >= v.nx())
                            j1 = v.nx() - 1;
                        return v(k1, i1, j1);
                    }
                    break;
                default:
                    break;
                }
            }
            break;
        
        default:
            break;
        }
    }

    template <int dim>
    void interp<dim>::cellLinearcorrect(int idim, double x, int &i1, int &i2, double &x1, double &x2){
        double midx = 0.5*( a[idim](i1) + a[idim](i2) );
        if (x < a[idim](0)){
            i1 = 0;
            i2 = 0;
            x1 = a[idim](0);
            x2 = 0.5*( a[idim](0) + a[idim](1) );
            return;
        }
        if (x > a[idim](a[idim].last())){
            if (idim == 0){
                i1 = v.nx() - 1;
                i2 = v.nx() - 1;
            }
            else if (idim == 1){
                i1 = v.ny() - 1;
                i2 = v.ny() - 1;
            }
            else if (idim == 2){
                i1 = v.nz() - 1;
                i2 = v.nz() - 1;
            }
            x1 = 0.5* ( a[idim](a[idim].last()-1) + a[idim](a[idim].last()) );
            x2 = a[idim](a[idim].last());
            return;
        }


        if (x < midx){
            x2 = midx;
            if (i1 != 0){
                x1 = 0.5*(a[idim](i1 - 1) + a[idim](i1));
                i2 = i1;
                i1 = i1 - 1;
            }
            else{
                x1 = a[idim](0);
                i2 = 0;
            }
        }
        else{
            x1 = midx;
            if (i2 != a[idim].last()){
                x2 = 0.5*(a[idim](i2) + a[idim](i2 + 1));
                //i1 = i2;
                //i2 = i2 + 1;
            }
            else{
                x2 = a[idim](a[idim].last());
                i2 = i1;
            }
        }
    }



    // /**
    //  * @brief interp1D is the class that it is used for 1D interpolation 
    //  * and for holding information for the axis for the higher dimension interpolations
    //  * 
    //  */
    // class interp1D{
    // public:
    //     //! This is a empty constructor 
    //     interp1D(){};
    //     /**
    //      * @brief Populates the positional vector. This is usefull when the space between
    //      * the interpolated is constant.
    //      * 
    //      * Starting from the Coordinate origin will add N times dx. 
    //      * For example setX(0,50, 11) will populate the vector x with the following numbers
    //      * 0 50 100 150 200 250 300 350 400 450 500
    //      * 
    //      * @param origin The starting coordinate. This corresponds to x[0]
    //      * @param dx the cell size
    //      * @param N The number of cells + 1. If the MODE is POINT then N is the number of ticks.
    //      * if MODE is CELL then N should be ncells + 1
    //      */
    //     void setX(double origin, double dx, int N);

    //     /**
    //      * @brief This simply sets the interpolation x vector equal to the input vector 
    //      * 
    //      * @param x_in 
    //      */
    //     void setX(std::vector<double> &x_in);
    //     /**
    //      * @brief This sets the values of the interpolation vector. 
    //      * It is the responsibility of the user to make sure that the inputs of X match the ones of the V and the methos type etc.
    //      * 
    //      * @param v_in 
    //      */
    //     void setV(std::vector<double> &v_in);
    //     void setModeMethod(MODE mode_in, METHOD method_in);
    //     double interpolate(double x_in);
    //     void findIIT(double x_in, int &i1, int &i2, double &t);
    //     void findIcell(double x_in, int &i);
    //     void clear();
    //     void dataFromFile(std::string filename);

    // private:
    //     std::vector<double> x;
    //     std::vector<double> v;
    //     MODE mode = MODE::INVALID_MODE;
    //     METHOD method = METHOD::INVALID_METHOD;
    //     int N;
    //     void findCell(int &i, int &ii, double x);
    //     bool isvalidated  = false;
    //     bool validateInputs();
    //     bool xConst = false;
    //     double dx = 0.0;
    //     double origin = 0.0;

    // };

    // void interp1D::clear(){
    //     x.clear();
    //     v.clear();
    //     mode = MODE::INVALID_MODE;
    //     method = METHOD::INVALID_METHOD;
    //     N = 0; 
    //     isvalidated = false;
    //     xConst = false;
    //     dx = 0.0;
    //     origin = 0.0;
    // }

    // void interp1D::setX(double origin_in, double dx_in, int n){
    //     isvalidated = false;
    //     x.clear();
    //     for(int i = 0; i < n; i++){
    //         x.push_back(origin_in + static_cast<double>(i)*dx_in);
    //     }
    //     xConst = true;
    //     origin = origin_in;
    //     dx = dx_in;
    // }
    // void interp1D::setX(std::vector<double> &x_in){
    //     isvalidated = false;
    //     xConst = false;
    //     x = x_in;
    // }

    // void interp1D::setV(std::vector<double> &v_in){
    //     isvalidated = false;
    //     v = v_in;
    // }

    // void interp1D::setModeMethod(MODE mode_in, METHOD method_in){
    //     isvalidated = false;
    //     mode = mode_in;
    //     method = method_in;
    // }

    // double interp1D::interpolate(double x_in){
    //     if (!isvalidated){
    //         isvalidated = validateInputs();
    //         if (!isvalidated)
    //             return std::numeric_limits<double>::quiet_NaN();
    //     }
    //     if (x_in < x[0])
    //         return v[0];
    //     else if (x_in > x[x.size()-1])
    //         return v[v.size()-1];
    //     else{
    //         int i1  = 0;
    //         int i2 = N;
    //         findCell(i1, i2, x_in);
    //         switch (mode)
    //         {
    //         case MODE::POINT:
    //             switch (method)
    //             {
    //             case METHOD::LINEAR:
    //                 {
    //                 double t = (x_in - x[i1])/(x[i2] - x[i1]);
    //                 return v[i1]*(1-t) + v[i2]*t;
    //                 }
    //                 break;
    //             case METHOD::NEAREST:
    //                 {
    //                 if (std::abs(x[i1] - x_in) <= std::abs(x[i2] - x_in))
    //                     return v[i1];
    //                 else
    //                     return v[i2];
    //                 }
    //                 break;
    //             default:
    //                 return std::numeric_limits<double>::quiet_NaN();
    //                 break;
    //             }
    //             break;

    //         case MODE::CELL:
    //             switch (method)
    //             {
    //             case METHOD::LINEAR:
    //                 {
    //                 double midx = 0.5*( x[i1] + x[i2] );
    //                 double x1, x2;
    //                 if (x_in < midx){
    //                     x2 = midx;
    //                     if (i1 != 0)
    //                         x1 = 0.5*(x[i1 - 1] + x[i1]);
    //                     else
    //                         return v[0];
                        
    //                     double t = (x_in - x1)/(x2 - x1);
    //                     return v[i1-1]*(1-t) + v[i1]*t;
    //                 }
    //                 else{
    //                     x1 = midx;
    //                     if (i2 != N)
    //                         x2 = 0.5*(x[i2] + x[i2 + 1]);
    //                     else
    //                         return v[v.size()-1];

    //                     double t = (x_in - x1)/(x2 - x1);
    //                     return v[i1]*(1-t) + v[i2]*t;
    //                 }
    //                 }
    //                 break;
    //             case METHOD::NEAREST:
    //                 return v[i1];
    //                 break;
    //             default:
    //                 return std::numeric_limits<double>::quiet_NaN();
    //                 break;
    //             }
    //             break;
    //         default:
    //             return std::numeric_limits<double>::quiet_NaN();
    //             break;
    //         }
    //     }
    // }

    // void interp1D::findCell(int &i, int &ii, double x_in){
    //     switch (xConst)
    //     {
    //     case true:
    //         {
    //             i = static_cast<int>((x_in - origin)/dx);
    //             ii = i + 1;
    //         }
    //         break;
    //     default:
    //         {
    //             int mid = (i + ii)/2;
    //             if (x_in <= x[mid])
    //                 ii = mid;
    //             else
    //                 i = mid;
                
    //             if (ii - i <= 1)
    //                 return;
    //             else
    //                 findCell(i, ii, x_in);
    //         }
    //         break;
    //     }
    // }

    // bool interp1D::validateInputs(){
    //     if (method == METHOD::INVALID_METHOD){
    //         std::cout << "The interpolation method is invalid" << std::endl;
    //         return false;
    //     }
    //     if (mode == MODE::INVALID_MODE){
    //         std::cout << "THe interpolation mode is invalid" << std::endl;
    //         return false;
    //     }
    //     if (mode == MODE::POINT && x.size() == v.size()){
    //         N = x.size()-1;
    //         return true;
    //     }
    //     if (mode == MODE::CELL && x.size() == v.size()+1){
    //         N = x.size()-1;
    //         return true;
    //     }
    //     std::cout << "The size of x is " << x.size() << " and the size of v is " <<  v.size() 
    //               << " which are not compatible with the selected MODE: " << mode <<  std::endl;
    //     return false;
    // }

    // void interp1D::findIIT(double x_in, int &i1, int &i2, double &t){
    //     if (x_in < x[0]){
    //         i1 = 0;
    //         i2 = 0;
    //         t = 0.0;
    //         return;
    //     }
    //     if (x_in > x[x.size()-1]){
    //         i1 = N;
    //         i2 = N;
    //         t = 1.0;
    //         return;
    //     }
    //     i1  = 0;
    //     i2 = N;
    //     findCell(i1, i2, x_in);
    //     t = (x_in - x[i1])/(x[i2] - x[i1]);
    // }

    // void interp1D::findIcell(double x_in, int &i){
    //     if (x_in < x[0]){
    //         i = 0;
    //         return;
    //     }
    //     if (x_in > x[x.size()-1]){
    //         i = N;
    //         return;
    //     }
    //     int i1  = 0;
    //     int i2 = N;
    //     findCell(i1, i2, x_in);
    //     i = i1;
    // }

    // void interp1D::dataFromFile(std::string filename){
    //     std::ifstream datafile(filename.c_str());
    //     if (!datafile.good()) {
    //         std::cout << "Can't open the file" << filename << std::endl;
    //     }
    //     else{
    //         MODE md;
    //         METHOD mthd;
    //         std::string line, tmp;
    //         bool isConstant;
    //         int has_values;
    //         getline(datafile, line);
    //         {
    //             std::istringstream inp(line.c_str());
    //             // Get mode
    //             inp >> tmp;
    //             if (tmp.compare("POINT") == 0)
    //                 md = MODE::POINT;
    //             else if (tmp.compare("CELL") == 0)
    //                 md = MODE::CELL;
    //             else
    //                 md = MODE::INVALID_MODE;

    //             //get method
    //             inp >> tmp;
    //             if (tmp.compare("LINEAR") == 0)
    //                 mthd = METHOD::LINEAR;
    //             else if (tmp.compare("NEAREST") == 0)
    //                 mthd = METHOD::NEAREST;
    //             else
    //                 mthd = METHOD::INVALID_METHOD; 

    //             inp >> tmp;
                
    //             if (tmp.compare("CONST") == 0)
    //                 isConstant = true;
    //             else
    //                 isConstant = false;

    //             inp >> has_values;
    //         }
            
    //         int N;
    //         double origin, dx;
    //         {
    //             getline(datafile, line);
    //             std::istringstream inp(line.c_str());
    //             if (isConstant){
    //                 inp >> origin;
    //                 inp >> dx;
    //                 inp >> N;
    //             }
    //             else{
    //                 inp >> N;
    //             }
    //         }
            
    //         std::vector<double> x,v;
    //         double vv, xx;
    //         if (!(has_values == 0 && isConstant)){
    //             for (int i = 0; i < N; i++){
    //                 getline(datafile, line);
    //                 std::istringstream inp(line.c_str());
    //                 if (isConstant){
    //                     inp >> vv;
    //                     v.push_back(vv); 
    //                 }
    //                 else{
    //                     inp >> xx;
    //                     x.push_back(xx);
    //                     if (has_values == 1){
    //                         inp >> vv;
    //                         v.push_back(vv);
    //                     }
    //                 }
    //             }
    //             if (md == MODE::CELL && !isConstant){
    //                 getline(datafile, line);
    //                 std::istringstream inp(line.c_str());
    //                 inp >> xx;
    //                 x.push_back(xx);
    //             }
    //         }
    //         if (isConstant){
    //             if (md == MODE::POINT)
    //                 setX(origin, dx, N);
    //             else if (md == MODE::CELL)
    //                 setX(origin, dx, N+1);
    //         }
    //         else{
    //             setX(x);
    //         }
    //         if (has_values == 1)
    //             setV(v);
    //         setModeMethod(md, mthd);
    //     }
    // }



/*
    class interp2D{
    public:
        interp2D(){}

        void setXY(std::vector<double> &x, std::vector<double> &y);
        void setXY(double xo, double yo, double dx, double dy, int nx, int ny);
        void setV(std::vector<std::vector<double>> &V_in);
        void setModeMethod(MODE mode_in, METHOD method_in);
        double interpolate(double x_in, double y_in);


    private:
        interp1D X;
        interp1D Y;
        std::vector<std::vector<double>> V;
        MODE mode = MODE::INVALID_MODE;
        METHOD method = METHOD::INVALID_METHOD;
        bool validateInputs();
        bool isvalidated = false;
    };

    void interp2D::setXY(std::vector<double> &x, std::vector<double> &y){
        X.setX(x);
        Y.setX(y);
    }

    void interp2D::setXY(double xo, double yo, double dx, double dy, int nx, int ny){
        X.setX(xo, dx, nx);
        Y.setX(yo, dy, ny);
    }

    void interp2D::setV(std::vector<std::vector<double>> &V_in){
        V =V_in;
    }

    void interp2D::setModeMethod(MODE mode_in, METHOD method_in){
        X.setModeMethod(mode_in, method_in);
        Y.setModeMethod(mode_in, method_in);
        mode = mode_in;
        method = method_in;
    }

    double interp2D::interpolate(double x_in, double y_in){
        if (!isvalidated){
            isvalidated = validateInputs();
            if (!isvalidated)
                return std::numeric_limits<double>::quiet_NaN();
        }
        int j1 = 0, j2 = 0, i1 = 0, i2 = 0;
        double tx = 0.0, ty = 0.0;
        X.findIIT(x_in, j1, j2, tx);
        Y.findIIT(y_in, i1, i2, ty);

        switch (mode)
        {
        case MODE::POINT:
            switch (method)
            {
            case METHOD::LINEAR:
                {
                    
                }
                break;
            case METHOD::NEAREST:
                {

                }
                break;
            default:
                return std::numeric_limits<double>::quiet_NaN();
                break;
            }
            break;
        case MODE::CELL:
            switch (method)
            {
            case METHOD::LINEAR:
                {
                    
                }
                break;
            case METHOD::NEAREST:
                {

                }
                break;
            default:
                return std::numeric_limits<double>::quiet_NaN();
                break;
            }
            break;
        default:
            return std::numeric_limits<double>::quiet_NaN();
            break;
        }
        return 0.0;
    }

    








    template<int dim>
    class interp{
    public:
        interp();

        void print();


    };


    template<int dim>
    interp<dim>::interp(){

    }

    template<int dim>
    void interp<dim>::print(){
        std::cout << "Interp Print" << std::endl;
    }
*/
    
    // void readData1D(std::string filename, interp1D &D){
    //     std::ifstream datafile(filename.c_str());
    //     if (!datafile.good()) {
    //         std::cout << "Can't open the file" << filename << std::endl;
    //     }
    //     else{
    //         MODE md;
    //         METHOD mthd;
    //         std::string line, tmp;
    //         bool isConstant;
    //         int has_values;
    //         getline(datafile, line);
    //         {
    //             std::istringstream inp(line.c_str());
    //             // Get mode
    //             inp >> tmp;
    //             if (tmp.compare("POINT") == 0)
    //                 md = MODE::POINT;
    //             else if (tmp.compare("CELL") == 0)
    //                 md = MODE::CELL;
    //             else
    //                 md = MODE::INVALID_MODE;

    //             //get method
    //             inp >> tmp;
    //             if (tmp.compare("LINEAR") == 0)
    //                 mthd = METHOD::LINEAR;
    //             else if (tmp.compare("NEAREST") == 0)
    //                 mthd = METHOD::NEAREST;
    //             else
    //                 mthd = METHOD::INVALID_METHOD; 

    //             inp >> tmp;
                
    //             if (tmp.compare("CONST") == 0)
    //                 isConstant = true;
    //             else
    //                 isConstant = false;

    //             inp >> has_values;
    //         }
            
    //         int N;
    //         double origin, dx;
    //         {
    //             getline(datafile, line);
    //             std::istringstream inp(line.c_str());
    //             if (isConstant){
    //                 inp >> origin;
    //                 inp >> dx;
    //                 inp >> N;
    //             }
    //             else{
    //                 inp >> N;
    //             }
    //         }
            
    //         std::vector<double> x,v;
    //         double vv, xx;
    //         if (!(has_values == 0 && isConstant)){
    //             for (int i = 0; i < N; i++){
    //                 getline(datafile, line);
    //                 std::istringstream inp(line.c_str());
    //                 if (isConstant){
    //                     inp >> vv;
    //                     v.push_back(vv); 
    //                 }
    //                 else{
    //                     inp >> xx;
    //                     x.push_back(xx);
    //                     if (has_values == 1){
    //                         inp >> vv;
    //                         v.push_back(vv);
    //                     }
    //                 }
    //             }
    //             if (md == MODE::CELL && !isConstant){
    //                 getline(datafile, line);
    //                 std::istringstream inp(line.c_str());
    //                 inp >> xx;
    //                 x.push_back(xx);
    //             }
    //         }
    //         if (isConstant){
    //             if (md == MODE::POINT)
    //                 D.setX(origin, dx, N);
    //             else if (md == MODE::CELL)
    //                 D.setX(origin, dx, N+1);
    //         }
    //         else{
    //             D.setX(x);
    //         }
    //         if (has_values == 1)
    //             D.setV(v);
    //         D.setModeMethod(md, mthd);
    //     }
    // }

    void writeCoordsValues(std::string filename, std::vector<std::vector<double>> &coords, std::vector<double> &v){
        std::ofstream out_file;
        out_file.open(filename.c_str());
        for (unsigned int i = 0; i < coords.size(); i++){
            for (unsigned j = 0; j < coords[i].size(); j++){
                out_file << coords[i][j] << " ";
            }
            out_file << v[i] << std::endl;
        }
        out_file.close();
    }
}

