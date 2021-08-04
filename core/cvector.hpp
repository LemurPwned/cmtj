#ifndef CVECTOR_H
#define CVECTOR_H
#include <random>
#include <stdio.h>
#include <vector>
template <typename T>
class CVector
{
public:
    T x, y, z;
    CVector()
    {
        this->x = 0.0;
        this->y = 0.0;
        this->z = 0.0;
    }
    CVector(std::vector<T> vec)
    {
        if (vec.size() != 3)
        {
            throw std::runtime_error("Failed to create vector -- passed list was not of len 3!");
        }
        this->x = vec[0];
        this->y = vec[1];
        this->z = vec[2];
    }

    CVector(T x, T y, T z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    CVector(const CVector &v)
    {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
    }

    CVector(std::normal_distribution<T> &distribution, std::default_random_engine &generator)
    {
        // the noise should be independent in each direction
        this->x = distribution(generator);
        this->y = distribution(generator);
        this->z = distribution(generator);
    }

    CVector &operator+=(const CVector &v)
    {
        this->x += v.x;
        this->y += v.y;
        this->z += v.z;
        return *this;
    }

    CVector &operator-=(const CVector &v)
    {
        this->x -= v.x;
        this->y -= v.y;
        this->z -= v.z;
        return *this;
    }
    CVector operator+(CVector v)
    {
        CVector res(
            x + v.x,
            y + v.y,
            z + v.z);

        return res;
    };

    CVector operator+(const CVector &v) const
    {
        CVector res(
            x + v.x,
            y + v.y,
            z + v.z);

        return res;
    };

    CVector operator+(const T &val) const
    {
        CVector res(
            x + val,
            y + val,
            z + val);
        return res;
    }
    CVector operator+(T &val) const
    {
        CVector res(
            x + val,
            y + val,
            z + val);
        return res;
    }

    CVector operator-(CVector v)
    {
        CVector res(
            x - v.x,
            y - v.y,
            z - v.z);

        return res;
    };
    CVector operator-(const CVector &v) const
    {
        CVector res(
            x - v.x,
            y - v.y,
            z - v.z);

        return res;
    };

    void operator=(CVector v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
    }
    bool operator==(CVector &v)
    {
        if (
            (x == v.x) && (y == v.y) && (y == v.z))
            return true;
        return false;
    };

    CVector operator*(T &val)
    {
        CVector res(
            x * val,
            y * val,
            z * val);
        return res;
    };

    CVector operator*(const T &val) const
    {
        const CVector res(
            x * val,
            y * val,
            z * val);
        return res;
    }

    CVector operator/(T val)
    {
        CVector res(
            x / val,
            y / val,
            z / val);
        return res;
    };
    T operator[](int &i)
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            return z;
    }

    T operator[](int &i) const
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            return z;
    }

    T operator[](const int &i) const
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            return z;
    }
    T length()
    {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }; // Magnitude
    void normalize()
    {
        T mag = this->length();
        if (mag != 0)
        {
            x = x / mag;
            y = y / mag;
            z = z / mag;
        }
    };
    void setX(T &vx)
    {
        this->x = vx;
    }
    void setY(T &vy)
    {
        this->y = vy;
    }
    void setZ(T &vz)
    {
        this->z = vz;
    }

    std::vector<T> tolist()
    {
        return {
            this->x, this->y, this->z};
    }

    operator std::string() const
    {
        return "[x:" + x + ", y:" + y + ", z:" + z + "]";
    }
};

#endif