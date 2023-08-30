//==================================================================================================
//  File:       matrix.hpp
//
//  Summary:    This header defines the matrix class template.
//==================================================================================================

#ifndef LINEAR_ALGEBRA_MATRIX_HPP
#define LINEAR_ALGEBRA_MATRIX_HPP

namespace LA_lib {

    template<typename T>
    class matrix {
    private:

        std::vector <T> buffer;
        size_t nrows;
        size_t ncols;

    public:

        using element_type = T;

// ==============================================================================
//                         MATRIX OBJECT INITIALIZATION
// ==============================================================================

        matrix() :
                buffer(0),
                nrows{0},
                ncols{0} {}

        template<typename N>
        matrix(N rows, N cols):
                buffer(rows * cols),
                nrows{static_cast<size_t>(rows)},
                ncols{static_cast<size_t>(cols)} {}

        template<typename N>
        matrix(std::initializer_list <T> &rhs, N rows, N cols):
                buffer{rhs},
                nrows{static_cast<size_t>(rows)},
                ncols{static_cast<size_t>(cols)} {}

        template<typename N>
        matrix(std::initializer_list <T> &&rhs, N rows, N cols):
                buffer{rhs},
                nrows{static_cast<size_t>(rows)},
                ncols{static_cast<size_t>(cols)} {}

        template<typename N>
        matrix(std::vector <T> &rhs, N rows, N cols):
                buffer{rhs},
                nrows{static_cast<size_t>(rows)},
                ncols{static_cast<size_t>(cols)} {}


// ==============================================================================
//                              MATRIX DATA ACCESS
// ==============================================================================

        T operator()(int ri, int ci) const
        {
            return buffer[(ri * ncols) + ci];
        }

        T &operator()(int ri, int ci)
        {
            return buffer[(ri * ncols) + ci];
        }

        T operator()(int i) const
        {
            assert(nrows == 1 || ncols == 1 && "2D object cannot be accessed by single index");

            return buffer[i];
        }

        T &operator()(int i)
        {
            assert(nrows == 1 || ncols == 1 && "2D object cannot be accessed by single index");

            return buffer[i];
        }

        T *data() { return buffer.data(); }


// ==============================================================================
//                            MATRIX PROPERTIES
// ==============================================================================

        [[nodiscard]] constexpr bool
        is_empty() const noexcept {
            if (nrows == 0 && ncols == 0)
                return true;
            else
                return false;
        }

        [[nodiscard]] constexpr size_t
        rows() const noexcept { return nrows; }

        [[nodiscard]] constexpr size_t
        columns() const noexcept { return ncols; }

        [[nodiscard]] constexpr size_t
        size() const noexcept { return (buffer.size()); }


// ==============================================================================
//                            MATRIX MODIFIERS
// ==============================================================================

        constexpr void
        swap_rows(int r1, int r2)
        {
            assert(r1 > 0 && r2 > 0 && "Invalid negative index");

            assert(r1 < nrows && r2 < nrows && "Row index is out of range");

            if (r1 != r2) {
                T aux;
                for (int i = 0; i < ncols; ++i) {
                    aux = this->operator()(r1, i);
                    this->operator()(r1, i) = this->operator()(r2, i);
                    this->operator()(r2, i) = aux;
                }
            }
        }

        constexpr void
        swap_columns(int c1, int c2)
        {
            assert(c1 > 0 && c2 > 0 && "Invalid negative index");

            assert(c1 < nrows && c2 < nrows && "Column index is out of range");

            if (c1 != c2) {
                T aux;
                for (int i = 0; i < nrows; ++i) {
                    aux = this->operator()(i, c1);
                    this->operator()(i, c1) = this->operator()(i, c2);
                    this->operator()(i, c2) = aux;
                }
            }
        }

        constexpr void
        resize(size_t rn, size_t cn)
        {
            assert(((rn == 0 && cn == 0) || (rn > 0 && cn > 0) && "Invalid matrix size"));

            if (rn != nrows || cn != ncols) {
                buffer.reserve(rn * cn);
                buffer.resize(rn * cn);
                buffer.shrink_to_fit();
                nrows = rn;
                ncols = cn;
            }
        }


// ==============================================================================
//                         ASSIGNMENT OPERATOR OVERLOADING
// ==============================================================================

        matrix<T> &operator=(matrix<T> const &rhs)
        {
            if (this != &rhs) {
                this->resize(rhs.nrows, rhs.ncols);
                buffer = rhs.buffer;
                nrows = rhs.nrows;
                ncols = rhs.ncols;
                return *this;
            } else {
                return *this;
            }
        }

        matrix<T> &operator=(std::initializer_list <T> const &rhs)
        {
            if (!std::empty(rhs)) {
                this->resize(1, rhs.size());
                buffer = rhs;
                nrows = 1;
                ncols = rhs.size();
                return *this;
            } else {
                this->resize(0, 0);
                buffer = rhs;
                nrows = 0;
                ncols = 0;
                return *this;
            }

        }

        matrix<T> &operator=(std::initializer_list <std::initializer_list<T>> const &rhs)
        {
            assert((check_equal_size_rows(rhs) && "The number of columns must be equal in every row"));

            if (is_list_empty(rhs)) {
                this->resize(0, 0);
                buffer = {};
                nrows = 0;
                ncols = 0;
                return *this;
            } else {
                this->resize(rhs.size(), rhs.begin()->size());

                int i = 0;
                for (auto const &v1: rhs) {
                    for (auto const &v2: v1) {
                        buffer[i] = v2;
                        ++i;
                    }
                }

                nrows = rhs.size();
                ncols = rhs.begin()->size();
                return *this;
            }
        }

        matrix<T> &operator=(matrix_view<T> const &rhs)
        {
            assert((nrows == rhs.rows() && ncols == rhs.columns()
                    && "Both matrices must have the same size"));

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = rhs(i,j);
                }
            }

            return *this;
        }

        template <class M>
        matrix<T> &operator=(matrix_transpose<M> const &rhs)
        requires
                contain_same_element_type<typename M::element_type, T>
        {
            assert((nrows == rhs.rows() && ncols == rhs.columns()
                    && "Both matrices must have the same size"));

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = rhs(i,j);
                }
            }

            return *this;
        }


// ==============================================================================
//                      ARITHMETIC OPERATORS OVERLOADING
// ==============================================================================

        // --------------------------- GENERAL ASSIGNMENT ---------------------------

        template<class M>
        matrix<T> &operator=(M const &rhs)
        {
            assert((nrows == rhs.rows() && ncols == rhs.columns()
                    && "Both matrices must have the same size"));

            for (int i = 0; i < nrows; ++i){
                for (int j = 0; j < ncols; ++j){
                    this->operator()(i,j) = rhs(i,j);
                }
            }

            return *this;
        }


        // --------------------------- ADDITION ---------------------------

        template <class M>
        matrix<T> & operator +=(const M &rhs)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            plus_equal_matrix<matrix<T>, M>(*this, rhs);
            return *this;
        }

        // --------------------------- SUBTRACTION ---------------------------

        template <class M>
        matrix<T> & operator -=(const M &rhs)
        requires
                contain_same_element_type<T, typename M::element_type>
        {
            minus_equal_matrix<matrix<T>, M>(*this, rhs);
            return *this;
        }

        matrix<T> operator-()
        {
            matrix<T> m{nrows, ncols};
            for (int i = 0; i < nrows; ++i) {
                for (int j = 0; j < ncols; ++j) {
                    m(i, j) = -this->operator()(i, j);
                }
            }
            return m;
        }


        // --------------------------- MATRIX PRODUCT ---------------------------

        template <class M>
        matrix<T> operator *=(const M & rhs)
        requires
        contain_same_element_type<typename M::element_type, T>
        {
            matrix<T> m;
            mult_matrix<matrix<T>, matrix<T>, M, T>(m, *this, rhs);
            return m;
        }


        // --------------------------- SCALAR PRODUCT ---------------------------

        matrix<T> & operator *=(const T &rhs)
        {
            matrix_times_equal_scalar<matrix<T>,T>(*this, rhs);
            return *this;
        }


        // --------------------------- DIVISION ---------------------------

        matrix<T> & operator /=(const T &rhs)
        {
            matrix_divide_equal_scalar<matrix<T>,T>(*this, rhs);
            return *this;
        }


        // --------------------------- COMPARISON ---------------------------

        matrix<T> operator==(matrix<T> const &rhs)
        {
            return (this->buffer == rhs.buffer);
        }

        matrix<T> operator!=(matrix<T> const &rhs)
        {
            return (this->buffer != rhs.buffer);
        }


// ==============================================================================
//                               SUBELEMENT ACCESS
// ==============================================================================

        matrix_view<T> submatrix(int ri, int rn, int ci, int cn)
        {
            return matrix_view(*this, ri, rn, ci, cn);
        }

        matrix_view<T> column(int i)
        {
            return matrix_view(*this, 0, nrows, i, 1);
        }

        matrix_view<T> row(int i)
        {
            return matrix_view(*this, i, 1, 0, ncols);
        }

        matrix_transpose<matrix<T>> t(){
            return matrix_transpose<matrix<T>>{*this};
        }
    };

}

#endif //LINEAR_ALGEBRA_MATRIX_HPP
