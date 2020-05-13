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
        void findIIT(double x_in, int &i1, int &i2, double &t)const;
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
        double operator()(int i)const;
        //! returns the index that access the last element of the tick values
        //! If ax is an axis object then ax(ax.last()) is equal to ax.x[ax.x.size()-1]
        int last()const;

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
        void findCell(int &i, int &ii, const double x)const;
    };

    void axis::reset(){
        x.clear();
        origin = 0.0;
        dx = 0.0;
        bConst = false;
        N = 0;
    }
    double axis::operator()(int i)const{
        return x[i];
    }
    int axis::last()const{
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

    void axis::findCell(int &i, int &ii, double x_in)const{
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

    void axis::findIIT(double x_in, int &i1, int &i2, double &t)const{
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


    /**
     * @brief This is a container class for the values.
     * 
     * When this object is constructed is empty. To start filling it with values
     * it has to be first resized using the #resize method.
     * After the class is resized it can be filled with values using the only set method
     */
    class GridValues{
    public:
        //! Empty constructor 
        GridValues(){};
        //! access operator in 3D 
        double operator()(int l, int r, int c)const;
        //! access operator in 2D
        double operator()(int r, int c)const;
        //! access operator in 1D 
        double operator()(int c)const;
        //! resize the data array. Unlike the usual usage of resize function in vector
        //! the resize will discard any existing data first 
        void resize(int nl, int nr, int nc);
        // sets the value of the l layer, r row and c column with the value v_in
        void set(int l, int r, int c, double v_in);
        //! Deletes the data of the object bringing it int the state to recieve otehr data
        void reset();
        //! returns the number of columns or values in the x direction
        int nx()const;
        //! returns the number of rows or values in the y direction
        int ny()const;
        //! returns the number of layers or values in the z direction
        int nz()const;
    private:
        std::vector<std::vector<std::vector<double>>> v;
    };

    double GridValues::operator()(int l, int r, int c)const{
        return v[l][r][c];
    }
    double GridValues::operator()(int r, int c)const{
        return v[0][r][c];
    }
    double GridValues::operator()(int c)const{
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

    int GridValues::nx()const{
        return v[0][0].size();
    }
    int GridValues::ny()const{
        return v[0].size();
    }
    int GridValues::nz()const{
        return v.size();
    }

    /**
     * @brief Main interpolation class.
     * 
     * This is a templated class where the dimention is the template parameter
     * 
     * @tparam dim is the dimension of the interpolation
     */
    template<int dim>
    class interp{
    public:
        // The constructor resets everything and allocates space for the axis.
        interp();
        //! This is the method to call to carry out the interpolation
        //! you have to pass as many coordinates are needed.
        double interpolate(double x, double y = 0, double z = 0)const;
        /**
         * @brief Set the Axis object
         * 
         * THis method is used when the ticks are not evenly spaced
         * 
         * @param idim is the index of the dimension. For the x, y and z this is 0,1 and 2 repsectively 
         * @param x_in the vector with the coordinates of the ticks
         */
        void setAxis(int idim, std::vector<double> &x_in);
        /**
         * @brief Set the Axis object when the ticks are evenly spaced
         * 
         * @param idim is the dimension index of the axis to set
         * @param origin_in is the starting point of the idim axis
         * @param dx_in is the space between the ticks
         * @param n the number of ticks in the idim axis
         */
        void setAxis(int idim, double origin_in, double dx_in, int n);
        //! sets the mode and the method
        void setModeMethod(MODE mode_in, METHOD method_in);

        /**
         * @brief Get the Data From File object
         * 
         * The data must be in the following format
         * 
         * 1st line:
         * MODE METHOD NX NY NZ
         * in case of 1 or 2D only NX or NX and NY are needed
         * 
         * 2nd line:
         * axisFileX axisFileY axisFileZ
         * This is a list of files that define the options for the 
         * X Y and Z axis.
         * 
         * For 1D repeat NX lines the values. 
         * 
         * For 2D repeat NY lines where each line containts NX values
         * 
         * For 3D repeat NZ times the NY x NX data.
         * 
         * Note that the number of values must be consistent with the number 
         * of axis values depending the selected MODE.
         * 
         * @param filename 
         */
        void getDataFromFile(std::string filename);
        //! Resets the class to the empty state.
        void reset();

    private:
        //! A vector of axis objects
        std::vector<axis> a;
        //! The container for the values that are used in the interpolation
        GridValues v;
        //! In 3D this holds the layer elevation information
        GridValues elev;
        //! This is the interpolation mode
        MODE mode = MODE::INVALID_MODE;
        //! THis holds the interpolation method
        METHOD method = METHOD::INVALID_METHOD;
        //! The interpolate method will call this when dim is 1
        double interp1D(double x)const;
        //! The interpolate method will call this when dim is 2
        double interp2D(double x, double y)const;
        //! The interpolate method will call this when dim is 2
        double interp3D(double x, double y, double z)const;
        //! A utility method which is used in MODE::CELL and METHOD::LINEAR
        //! Finds the correct cells to interpolate the value for the given input
        void cellLinearcorrect(int idim, double x, int &i1, int &i2, double &x1, double &x2)const;
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
        elev.reset();
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
    double interp<dim>::interpolate(double x, double y, double z)const{
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
    double interp<dim>::interp1D(double x)const{
        double outcome;
        int i1, i2;
        double t;
        a[0].findIIT(x,i1, i2, t);
        
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
    double interp<dim>::interp2D(double x, double y)const{
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
    double interp<dim>::interp3D(double x, double y, double z)const{
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
    void interp<dim>::cellLinearcorrect(int idim, double x, int &i1, int &i2, double &x1, double &x2)const{
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

    //! This is a utility function which is used to print the results of the interpolation to a file so that it can be read in R forexample
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

