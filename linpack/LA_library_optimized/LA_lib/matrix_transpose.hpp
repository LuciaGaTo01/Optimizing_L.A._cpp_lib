//==================================================================================================
//  File:       matrix_transpose.hpp
//
//  Summary:    This header defines the matrix_transpose class template.
//==================================================================================================

#ifndef LINEAR_ALGEBRA_MATRIX_TRANSPOSE_HPP
#define LINEAR_ALGEBRA_MATRIX_TRANSPOSE_HPP

namespace LA_lib{

    template <class M>
    class matrix_transpose
    {
    private:
        M matt;
        size_t nrows;
        size_t ncols;

    public:

        using T = typename M::element_type;
        using element_type = T;

// ==============================================================================
//                          MATRIX TRANSPOSE INITIALIZATION
// ==============================================================================

        explicit matrix_transpose(matrix<T> &m):
                matt{m},
                nrows{m.columns()},
                ncols{m.rows()}
        {}

        explicit matrix_transpose(matrix_view<T> &m):
                matt{m},
                nrows{m.columns()},
                ncols{m.rows()}
        {}


// ==============================================================================
//                            MATRIX VIEW DATA ACCESS
// ==============================================================================

        T operator()(int ri, int ci) const
        {
            if (ri >= nrows || ci >= ncols)
                throw std::out_of_range("Index is out of range");

            return matt(ci, ri);
        }

        T &operator()(int ri, int ci)
        {
            if (ri >= nrows || ci >= ncols)
                throw std::out_of_range("Index is out of range");

            return matt(ci, ri);
        }

        T operator()(int i) const { return matt(i); }

        T &operator()(int i) { return matt(i); }


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
        size() const noexcept { return (nrows*ncols); }


// ==============================================================================
//                            MATRIX MODIFIERS
// ==============================================================================

        constexpr void
        swap_rows(int r1, int r2)
        {
            matt.swap_columns(r1, r2);
        }

        constexpr void
        swap_columns(int c1, int c2)
        {
            matt.swap_rows(c1,c2);
        }

        matrix_transpose<M> operator-()
        {
            matrix_transpose<M> m{-matt};
            return m;
        }


// ==============================================================================
//                      ARITHMETIC OPERATORS OVERLOADING
// ==============================================================================

        // --------------------------- COMPARISON ---------------------------

        matrix_transpose<M> operator==(matrix_transpose<M> const &rhs)

        {
            return matt == rhs.matt;
        }

        matrix_transpose<M> operator!=(matrix_transpose<M> const &rhs)
        {
            return matt != rhs.matt;
        }


    };



}


#endif //LINEAR_ALGEBRA_MATRIX_TRANSPOSE_HPP
