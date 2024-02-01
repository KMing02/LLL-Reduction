#include <iostream>
#include <vector>
#include <sstream>
#include <numeric>
#include <string>

// Function to display a vector of vectors
void displayVectors(const std::vector<std::vector<double> >& vectors) {
    for (const auto& vector : vectors) {
        std::cout << "Vector:";
        for (const auto& element : vector) {
            std::cout << " " << element;
        }
        std::cout << std::endl;
    }
}

// Function to display a vector
void displayVector(const std::vector<double>& vector) {
    for (const auto& element : vector) {
        std::cout << " " << element;
    }
}

// Define dot product
double dotProduct(const std::vector<double>& v1, const std::vector<double>& v2) {
    double product = std::inner_product(v1.begin(),v1.end(),v2.begin(),0.0);
    return product;
}

// Define vector subtraction, the returned vector is the result of elementwise u - v
std::vector<double> vectorSubtraction(std::vector<double>& u, std::vector<double>& v) {
    std::vector<double> resultVector(u.size());
    for (std::size_t i = 0; i < u.size(); ++i) {
        resultVector[i] = u[i] - v[i];
    }
    return resultVector;
}

// Define vector subtraction, the returned vector is the result of elementwise u - v
std::vector<double> vectorAddition(std::vector<double>& u, std::vector<double>& v) {
    std::vector<double> resultVector(u.size());
    for (std::size_t i = 0; i < u.size(); ++i) {
        resultVector[i] = u[i] + v[i];
    }
    return resultVector;
}

// Define vector subtraction, the returned vector is the result of elementwise u - v
std::vector<double> scalarMultiplication(std::vector<double>& u, double scalar) {
    std::vector<double> resultVector(u.size());
    for (std::size_t i = 0; i < u.size(); ++i) {
        resultVector[i] = u[i] * scalar;
    }
    return resultVector;
}

std::vector<double> vectorProjection(std::vector<double>& u, std::vector<double>& v) {
    std::vector<double> resultVector(u.size());
    double scalar = dotProduct(u,v) / dotProduct(u,u);
    resultVector =  scalarMultiplication(u, scalar);
    return resultVector;
}

// Function to perform Gram-Schmidt orthogonalization without normalization
std::vector<std::vector<double> > gramSchmidt(std::vector<std::vector<double> >& vectors) {
    std::vector<std::vector<double> > GSvectors(vectors.size());
    for (std::size_t i = 0; i < vectors.size(); ++i) {
        GSvectors[i] = vectors[i];
        for (std::size_t j = 0; j < i; ++j) {
            std::vector<double> projectionVector = vectorProjection(GSvectors[j], vectors[i]);
            GSvectors[i] = vectorSubtraction(GSvectors[i], projectionVector);
        }
    }
    return GSvectors;
}

/*  Lenstra–Lenstra–Lovász lattice basis reduction algorithm
    vectors:
    GSvectors:
    delta:
*/
void LLL(std::vector<std::vector<double> >& vectors, double delta) {
    std::vector<std::vector<double> > gs = gramSchmidt(vectors);

    int n = vectors.size();
    gs = gramSchmidt(vectors);
    std::vector<std::vector<double> > mu(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mu[i][j] = dotProduct(vectors[i], gs[j]) / dotProduct(gs[j], gs[j]);
        }
    }
    int k = 1;

    std::cout << "the gs vectors are" << std::endl;
    displayVectors(gs);
    std::cout << mu[1][0] << std::endl;

    while (k < n) {
        for (int j = k-1; j >= 0; j--) {
            std::vector<double> temp = scalarMultiplication(vectors[j], std::round(mu[k][j]));
            std::cout << k << j << mu[k][j] << std::endl;
            displayVector(temp);

            vectors[k] = vectorSubtraction(vectors[k], temp);
            std::cout << "updated vectors" << std::endl;
            displayVectors(vectors);
            gs = gramSchmidt(vectors);
            mu[k][j] = dotProduct(vectors[k], gs[j]) / dotProduct(gs[j], gs[j]);
        }
        if (dotProduct(gs[k], gs[k]) > (delta - mu[k][k-1] * mu[k][k-1]) * (dotProduct(gs[k-1], gs[k-1]))) {
            k++;
        } else {
            std::cout << "vectors before swap" << std::endl;

            std::vector<double> temp = vectors[k];
            vectors[k] = vectors[k-1];
            vectors[k-1] = temp;
            gs = gramSchmidt(vectors);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    mu[i][j] = dotProduct(vectors[i], gs[j]) / dotProduct(gs[j], gs[j]);
                }
            }
            k = std::max(1,k-1);
            std::cout << "vectors after swap" << std::endl;
            displayVectors(vectors);
            std::cout << "gsvectors after swap" << std::endl;
            displayVectors(gs);
        }
    }
}

double euclideanDistance(std::vector<double> & vector) {
    int dim = vector.size();
    double dist = 0.0;
    for (int i = 0; i < dim; i++) {
        dist += vector[i] * vector[i];
    }
    dist /= dim;
    return dist;
}

double findShortestVector(std::vector<std::vector<double> >& vectors) {
    int dim = vectors.size();
    double minlength = euclideanDistance(vectors[0]);
    for (int i = 1; i < dim; i++) {
        if (euclideanDistance(vectors[i]) < minlength) {
            minlength = euclideanDistance(vectors[i]);
        }
    }
    return minlength;
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./runme [vector1] [vector2] ... [vectorN]" << std::endl;
        return 1;
    }

    // Initialize a vector to store vectors
    std::vector<std::vector<double> > vectors;

    // Initialize a vector to store Gram-Schmidt vectors
    std::vector<std::vector<double> > GSvectors(vectors.size());

    // Loop through command-line arguments (skip the program name at index 0)
    for (int i = 1; i < argc; ++i) {
        // Create a stream to parse the vector string
        std::istringstream vectorStream(argv[i]);

        // Parse the vector elements and store them in a vector
        std::vector<double> vectorElements;
        double element;
        while (vectorStream >> element) {
            vectorElements.push_back(element);
        }

        // Add the vector to the vector of vectors
        vectors.push_back(vectorElements);
    }

    // Display the vectors
    std::cout << "initial basis" << std::endl;
    displayVectors(vectors);

    std::vector<std::vector<double> > testGS(vectors.size());
    testGS = gramSchmidt(vectors);

    displayVectors(testGS);

    LLL(vectors, 0.75);
    displayVectors(vectors);

    // Specify the file name
    std::string filename = "example.txt";

    // Create an input stream object
    std::ifstream inputFile(filename);

    // We approach the answer by firstly using LLL algorithm for reducing the basis vectors
    // and then pick the shortest vector from the list for reduced basis vectors.
    
    return 0;
}
