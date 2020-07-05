#ifndef CVECTOR_H
#define CVECTOR_H 
#include <stdio.h>
#include <random>
class CVector
{
public:
    double x, y, z;          // Order list of 4 elements |x|y|z|w|
    static char sBuffer[38]; // holds the string of given vector
    char *toString();        // Return sBuffer with (x,y,z,w) values
    CVector()
    {
        this->x = 0.0;
        this->y = 0.0;
        this->z = 0.0;
    } // zero Vector Constructor
    CVector(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    } // Constructor
    CVector(const CVector &v)
    {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
    } // Copy Vector Constructor

    CVector(std::normal_distribution<double> &distribution, std::default_random_engine &generator)
    {
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
    }; // Destructor
        // Addition
    CVector operator-(CVector v)
    {
        CVector res(
            x - v.x,
            y - v.y,
            z - v.z);

        return res;
    };
    double dotProduct(CVector &);

    void operator=(CVector v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
    } // Copy Vector
    bool operator==(CVector &v)
    {
        if (
            (x == v.x) && (y == v.y) && (y == v.z))
            return true;
        return false;
    };

    CVector operator*(double val)
    {
        CVector res(
            x * val,
            y * val,
            z * val);
        return res;
    };

    CVector operator/(double val)
    {
        CVector res(
            x / val,
            y / val,
            z / val);
        return res;
    };
    double operator[](int i)
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            return z;
    } // Scalar Multiplication
    double length()
    {
        return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    }; // Magnitude
    void normalize()
    {
        double mag = this->length();
        x = x / mag;
        y = y / mag;
        z = z / mag;
    };
    void setX(double vx)
    {
        this->x = vx;
    }
    void setY(double vy)
    {
        this->y = vy;
    }
    void setZ(double vz)
    {
        this->z = vz;
    }
};

#endif