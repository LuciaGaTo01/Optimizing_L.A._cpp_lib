

#ifndef ALINEAR_ALGEBRA_DEFINITIONS_HPP_DEFINED
#define ALINEAR_ALGEBRA_DEFINITIONS_HPP_DEFINED

namespace LA_lib{

// template <typename T>
// class submat{
//     private:

//     matrix<T> &m;
//     int ri;
//     int rn;
//     int ci;
//     int cn;
//     int n_elems;

//     public:

//     submat():
//         ri {0},
//         rn {0},
//         ci {0},
//         cn {0},
//         n_elems {0}
//     {}

//     submat(matrix<T> &ma, int icol, int ncol, int irow, int nrow):
//         m {ma},
//         ri{irow},
//         rn{nrow},
//         ci{icol},
//         cn{ncol},
//         n_elems{nrow*ncol}
//     {}

//     constexpr int
//     rows() const noexcept
//     {
//         return this->rn;
//     }

//     constexpr int
//     columns() const noexcept
//     {
//         return this->rc;
//     }

//     constexpr int
//     size() const noexcept
//     {
//         return this->n_elems;
//     }

//     constexpr auto
//     operator () (int r, int c)
//     {
//         return m(ri+r,c+ci);
//     }

//     template <typename TT>
//     constexpr submat&
//     operator =(submat<TT> const& subminit)
//     {
//         for (int i = 0; i < rn; ++i){
//             for (int j = 0; j < cn; ++j){
//                 m(i+ri,j+ci) = subminit(i,j);
//             }
//         }

//         return *this;
//     }

//     template <typename TT>
//     constexpr submat
//     operator +(submat<TT> const& m1, submat<TT> const& m2)
//     {
//         for (int i = 0; i < rn; ++i){
//             for (int j = 0; j < cn; ++j){
//                 m(i+ri,j+ci) = m1(i,j) + m2(i,j);
//             }
//         }

//         return *this;
//     }

//     constexpr submat
//     operator -()
//     {
//         for (int i = 0; i < rn; ++i){
//             for (int j = 0; j < cn; ++j){
//                 m(i+ri,j+ci) = - m(i+ri,j+ci);
//             }
//         }

//         return *this;
//     }

//     template <typename TT>
//     constexpr submat
//     operator -(submat<TT> const& m1,submat<TT> const& m2)
//     {
//         this->resize(m1.rows(), m1.columns());

//         for (int i = 0; i < n_rows; ++i){
//             for (int j = 0; j < n_cols; ++j){
//                 m(ri+i,ci+j) = m1(i,j) - m2(i,j);
//             }
//         }

//         return *this;
//     }

//     template <typename TT>
//     constexpr submat
//     operator *(submat<TT> const& m1,submat<TT> const& m2)
//     {
//         this->resize(m1.rows(), m2.columns());

//         for (int i = 0; i < n_rows; ++i){
//             for (int j = 0; j < n_cols; ++j){

//                 TT aux = 0;
//                 for (int k = 0; k < m1.rows(); ++k){
//                     aux += m1(i,k) * m2(k,j);
//                 }
//                 m(ri+i,ci+j) = aux;
//             }
//         }

//         return *this;
//     }

//     template <typename TT>
//     constexpr submat
//     operator *(submat<TT> const& m1, TT const& val)
//     {
//         this->resize(m1.rows(), m1.columns());

//         for (int i = 0; i < n_rows; ++i){
//             for (int j = 0; j < n_cols; ++j){
//                 m(ri+i,ci+j) = m1(i,j) * val;
//             }
//         }

//         return *this;
//     }

//     template <typename TT>
//     constexpr submat
//     operator *(TT const& val, submat<TT> const& m1)
//     {
//         this->resize(m1.rows(), m1.columns());

//         for (int i = 0; i < n_rows; ++i){
//             for (int j = 0; j < n_cols; ++j){
//                 m(ri+i,ci+j) = val * m1(i,j);
//             }
//         } 

//         return *this;
//     }

//     template <typename TT>
//     constexpr submat
//     operator /(submat<TT> const& m1, TT const& val)
//     {
//         this->resize(m1.rows(), m1.columns());

//         for (int i = 0; i < n_rows; ++i){
//             for (int j = 0; j < n_cols; ++j){
//                 m(ri+i,ci+j) = m1(i,j)/val;
//             }
//         }

//         return *this;
//     }


// };



template <typename T>
class matrix{
    private:

    using matrix_type = std::vector<T>;

    matrix_type m_elems;
    int n_rows;
    int n_cols;
    int n_elems;

    public:

    matrix():
        m_elems(0),
        n_rows{0},
        n_cols{0},
        n_elems{0}
    {}

    matrix(T rows, T cols): 
        m_elems(rows*cols), 
        n_rows{rows},
        n_cols{cols},
        n_elems{rows*cols}
    {}

    template <class U>
    constexpr matrix&
    operator =(std::initializer_list<U> linit)
    {
        m_elems = linit;
        return *this;
    }

    template <typename TT>
    constexpr matrix&
    operator =(matrix<TT> const& minit)
    {
        m_elems = minit.m_elems;
        n_rows = minit.n_rows;
        n_cols = minit.n_cols;
        n_elems = minit.n_elems;
        return *this;
    }

    constexpr int
    rows() const noexcept
    {
        return this->n_rows;
    }

    constexpr int
    columns() const noexcept
    {
        return this->n_cols;
    }

    constexpr int
    size() const noexcept
    {
        return this->n_elems;
    }

    constexpr void
    resize(int rn, int cn)
    {
        m_elems.reserve(rn*cn);
        m_elems.resize(rn*cn);
        n_rows = rn;
        n_cols = cn;
        n_elems = rn*cn;
    }

    constexpr void
    resize_rows(int rn)
    {
        m_elems.reserve(rn*n_cols);
        m_elems.resize(rn*n_cols);
        n_rows = rn;
        n_elems = rn*n_cols;
    }

    constexpr void
    resize_columns(int cn)
    {
        m_elems.reserve(n_rows*cn);
        m_elems.resize(n_rows*cn);
        n_cols = cn;
        n_elems = n_rows*cn;
    }

    constexpr auto
    operator () (int ri, int ci)
    {
        return m_elems[(ri*n_cols) + ci];
    }

    constexpr auto
    operator () (int i)
    {
        return m_elems[i];
    }

    constexpr matrix
    operator -()
    {
        for (int i = 0; i < n_rows; ++i){
            for (int j = 0; j < n_cols; ++j){
                m_elems[(i*n_cols)+j] = - m_elems[(i*n_cols)+j];
            }
        }

        return *this;
    }

    // template <typename TT>
    // constexpr submat
    // submatrix(int ri, int rn, int ci, int cn)
    // {
    //     return submat(*this, ci, cn, ri, rn);
    // }

};

template <typename TT>
constexpr matrix<TT>
operator +(matrix<TT> const& m1, matrix<TT> const& m2)
{
    matrix<TT> m(m1.rows(), m1.columns());

    for (int i = 0; i < m1.rows(); ++i){
        for (int j = 0; j < m1.columns(); ++j){
            m(i,j) = m1(i,j) + m2(i,j);
        }
    }

    return m;

}

template <typename TT>
constexpr matrix<TT>
operator -(matrix<TT> const& m1,matrix<TT> const& m2)
{
    matrix<TT> m(m1.rows(), m1.columns());

    for (int i = 0; i < m1.rows(); ++i){
        for (int j = 0; j < m1.columns(); ++j){
            m(i,j) = m1(i,j) - m2(i,j);
        }
    }

    return m;
}

template <typename TT>
constexpr matrix<TT>
operator *(matrix<TT> const& m1,matrix<TT> const& m2)
{
    matrix<TT> m(m1.rows(), m2.columns());

    for (int i = 0; i < m1.rows(); ++i){
        for (int j = 0; j < m2.columns(); ++j){

            TT aux = 0;
            for (int k = 0; k < m1.rows(); ++k){
                aux += m1(i,k) * m2(k,j);
            }
            m(i,j) = aux;
        }
    }

    return m;
}

template <typename TT>
constexpr matrix<TT>
operator *(matrix<TT> const& m1, TT const& val)
{
    matrix<TT> m(m1.rows(), m1.columns());

    for (int i = 0; i < m1.rows(); ++i){
        for (int j = 0; j < m1.columns(); ++j){
            m(i,j) = m1(i,j) * val;
        }
    }

    return m;
}

template <typename TT>
constexpr matrix<TT>
operator *(TT const& val, matrix<TT> const& m1)
{
    matrix<TT> m(m1.rows(), m1.columns());

    for (int i = 0; i < m1.rows(); ++i){
        for (int j = 0; j < m1.columns(); ++j){
            m(i,j) = val * m1(i,j);
        }
    } 

    return m;
}

template <typename TT>
constexpr matrix<TT>
operator /(matrix<TT> const& m1, TT const& val)
{
    matrix<TT> m(m1.rows(), m1.columns());

    for (int i = 0; i < m1.rows(); ++i){
        for (int j = 0; j < m1.columns(); ++j){
            m(i,j) = m1(i,j)/val;
        }
    }

    return m;
}

}


#endif