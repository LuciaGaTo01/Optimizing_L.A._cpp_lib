
#ifndef ALINEAR_ALGEBRA_DEFINITIONS_HPP_DEFINED
#define ALINEAR_ALGEBRA_DEFINITIONS_HPP_DEFINED

namespace LA_lib{

template <typename T>
class matrix_view;


template <typename T>
class matrix
{
    private:

        std::vector<T> buffer;
        int nrows;
        int ncols;

        template <typename U>
        friend class matrix_view;

    public:

        matrix():
            buffer(0),
            nrows{0},
            ncols{0}
        {}

        matrix(T rows, T cols):
            buffer(rows*cols),
            nrows{rows},
            ncols{cols}
        {}

        template <typename U>
        constexpr matrix&
        operator =(matrix<U> const& rhs)
        {
            buffer = rhs.buffer;
            nrows = rhs.nrows;
            ncols = rhs.ncols;
            return *this;
        }

        constexpr int
        rows() const noexcept{ return nrows; }

        constexpr int
        columns() const noexcept{ return ncols; }

        constexpr int
        size() const noexcept{ return (nrows*ncols); }

        constexpr void
        resize(int rn, int cn)
        {
            buffer.reserve(rn*cn);
            buffer.resize(rn*cn);
            nrows = rn;
            ncols = cn;
        }

        T operator () (int i) const{ return buffer[i]; }

        T &operator () (int i){ return buffer[i]; } 

        T operator () (int ri, int ci) const{ return buffer[(ri*ncols) + ci]; }

        T &operator () (int ri, int ci){ return buffer[(ri*ncols) + ci]; } 

        matrix_view<T> submatrix(int ri, int rn, int ci, int cn)
        {
            return matrix_view(*this, ri, rn, ci, cn);
        }

        T * data(){
            return buffer.data();
        }

};


template <typename T>
class matrix_view{
    private:

        matrix<T> &mview;
        int row_index;
        int row_size;
        int column_index;
        int column_size;

    public:

        matrix_view(matrix<T> &m, int ri, int rn, int ci, int cn): 
            mview{m}, 
            row_index{ri}, row_size{rn}, 
            column_index{ci}, column_size{cn}
        {}

        constexpr int
        rows() const noexcept{ return row_size; }

        constexpr int
        columns() const noexcept{ return column_size; }

        constexpr int
        size() const noexcept{ return (row_size*column_size); }

        T operator () (int i) const{ return mview(i+row_index); }

        T &operator () (int i){ return mview(i+row_index); } 

        T operator () (int ri, int ci) const{ return mview(ri+row_index, ci+column_index); }

        T &operator () (int ri, int ci){ return mview(ri+row_index, ci+column_index); } 
        
        matrix_view<T> & operator +=(const matrix_view<T> &rhs)
        {
            for (int i = 0; i < row_size; ++i){
                for (int j = 0; j < column_size; ++j){
                    this->operator()(i,j) += rhs(i,j);
                }
            }

            return *this;
        }

        matrix_view<T> & operator *=(const T &rhs)
        {
            for (int i = 0; i < row_size; ++i){
                for (int j = 0; j < column_size; ++j){
                    this->operator()(i,j) *= rhs;
                }
            }

            return *this;
        }

        void swap_columns(int ci, int cj)
        {
            T aux = 0;
            if (ci != cj){
                for (int i = 0; i < row_size; ++i){
                    aux = this->operator()(i,ci);
                    this->operator()(i,ci) = this->operator()(i,cj);
                    this->operator()(i,cj) = aux;
                }
            }

        }

        friend matrix<T> operator *(const matrix_view &lhs, const matrix_view &rhs){
            matrix<T> mat(lhs.rows(), rhs.columns());

            T aux;
            for (int i = 0; i < lhs.rows(); ++i){
                for (int j = 0; j < rhs.columns(); ++j){
                    aux = 0;
                    for (int k = 0; k < lhs.columns(); ++k){
                        aux += lhs(i,k) * rhs(k, j);
                    }
                    mat(i,j) = aux;
                }
            }

            return mat;

        }

};

template <class T>
void PRINT (T &mat){
    std::cout << "\n\nSize: " << mat.rows() << "x" << mat.columns() << "\n";
    std::cout << "---------------------------------------------------------\n";
    for (int i = 0; i < mat.rows(); ++i){
        for (int j = 0; j < mat.columns(); ++j){
            std::cout << mat(i,j) << "\t";
        }
        std::cout << "\n";
    }
}

}


#endif