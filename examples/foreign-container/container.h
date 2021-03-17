#ifndef LINEARALGEBRA_CONTAINER_H
#define LINEARALGEBRA_CONTAINER_H


template <typename T=double>
class container {
    std::vector<T> m_data;
public:

    container(size_t n) : m_data(n) {}

    using value_type = typename decltype(m_data)::value_type;

    std::vector<value_type> &data() { return m_data; }

    const std::vector<value_type> &data() const { return m_data; }

    void fill(value_type a) { std::for_each(m_data.begin(), m_data.end(), [a](value_type &y) { y = a; }); }

    void scal(value_type a) { std::for_each(m_data.begin(), m_data.end(), [a](value_type &y) { y *= a; }); }

    void axpy(value_type a, const container &x) {
        std::transform(x.m_data.begin(), x.m_data.end(),
                       m_data.begin(),
                       m_data.begin(),
                       [a](const value_type &xx, const value_type &yy) { return yy + a * xx; });
    }

    value_type dot(const container &x) const {
        return std::inner_product(m_data.begin(), m_data.end(), x.m_data.begin(), static_cast<value_type>(0));
    }

    std::map<size_t, value_type> select_max_dot(size_t n, const container& y) const {throw std::logic_error("unimplemented");}
};


#endif //LINEARALGEBRA_CONTAINER_H
