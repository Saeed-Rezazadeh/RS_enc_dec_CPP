#include <iostream>
#include <cmath>
#include <vector>
#include "GaloisField.h"
#include "GaloisFieldElement.h"
#include "GaloisFieldPolynomial.h"

using namespace std;
using namespace galois;

GaloisFieldElement determinant(GaloisFieldElement **, int n);

unsigned int prim_poly[] = {1, 0, 1, 1, 0, 1, 0, 0, 1};  // Primitive polynomial p(x) coefficients for GF(2 ^ 8)
//unsigned int prim_poly [] = {1 , 1 , 0 , 0 , 1};  // Primitive polynomial p(x) coefficients for GF(2 ^ 4)
unsigned int message[] = {116, 115, 101, 84};

int m = 8; // where m determines the order of the field.
/*
  A Galois Field of type GF(2^m)
*/

GaloisField galois_field(m, prim_poly);


int main(int argc, char *argv[]) {
    int n = 12; // The total size of the block message (k+redundant)
    int k = 4; //   k is the size of the message

    // Find the message polynomial such that P(X) = X^{2*t} * a(X) where a(x) is polynomial of degree k - 1
    // whose coefficients are the message to be encoded.
    GaloisFieldElement messageElements[n];
    int messageCounter = 0;
    for (int i = 0; i < n; i++) {
        if (i < n - k) {
            messageElements[i] = GaloisFieldElement(&galois_field, 0);
        } else {
            messageElements[i] = GaloisFieldElement(&galois_field, message[messageCounter++]);
        }
    }// for statement

    // Find the message polynomial
    GaloisFieldPolynomial messagePolynomial = GaloisFieldPolynomial(&galois_field, n - 1, messageElements);

    // Find the generator polynomial such that all elements of GF(2^m) are its roots
    GaloisFieldElement generatorElement[1] = {GaloisFieldElement(&galois_field, 1)};
    GaloisFieldPolynomial generatorPolynomial = GaloisFieldPolynomial(&galois_field, 0, generatorElement); // g_0(X) = 1

    GaloisFieldElement intermediateElements[2] = {GaloisFieldElement(&galois_field, galois_field.alpha(1)),
                                                  GaloisFieldElement(&galois_field, 1)};
    GaloisFieldPolynomial intermediatePolynomial;
    // Find g(X) = (X + alpha) * ... * (X + alpha ^ {n - k - 1})
    for (int i = 0; i < n - k; i++) {
        intermediateElements[0] = GaloisFieldElement(&galois_field, galois_field.alpha(i + 1));
        generatorPolynomial *= GaloisFieldPolynomial(&galois_field, 1, intermediateElements);
    } // for statement
    GaloisFieldPolynomial remainderPolynomial = messagePolynomial % generatorPolynomial;

    // Find the encodedPolynomial
    GaloisFieldPolynomial encodedPolynomial = messagePolynomial + remainderPolynomial;

    GaloisFieldPolynomial receivedPolynomial = encodedPolynomial;
    // Let's assume the following bits have been modified due to channel noise
    receivedPolynomial[4] = 69;
    receivedPolynomial[5] = 38;
    receivedPolynomial[9] = 121;
    receivedPolynomial[11] = 170;

    // Find the syndromes by evaluating the receivedPolynomial at alpha ^ i for i = 1 , ... , 2t
    GaloisFieldElement syndromeElements[n - k];
    for (int i = 0; i < n - k; i++) {
        syndromeElements[i] = GaloisFieldElement(&galois_field, 0);
    }
    for (int i = 1; i <= n - k; i++) {
        for (int j = 0; j < n; j++) {
            syndromeElements[n - k - i] =
                    syndromeElements[n - k - i] + galois_field.alpha(j * i) * receivedPolynomial[j];

        }// for statement j
    }// for statement i
    GaloisFieldPolynomial syndromePolynomial = GaloisFieldPolynomial(&galois_field, n - k - 1, syndromeElements);
    // The decoder is based on Euclidean Algorithm described in
    // Error Control Coding: Fundamentals and Applications by Lin and Constello 2nd Ed.

    GaloisFieldElement unitElement[1] = {GaloisFieldElement(&galois_field, 1)};
    GaloisFieldPolynomial Zi_2 = GaloisFieldPolynomial(&galois_field, 0, unitElement);
    Zi_2 <<= n - k;

    GaloisFieldPolynomial Zi_1 = syndromePolynomial;

    GaloisFieldElement zeroElement[1] = {GaloisFieldElement(&galois_field, 0)};
    GaloisFieldPolynomial sigma_i_2 = GaloisFieldPolynomial(&galois_field, 0, zeroElement);

    GaloisFieldPolynomial sigma_i_1 = GaloisFieldPolynomial(&galois_field, 0, unitElement);


    GaloisFieldPolynomial Zi;
    GaloisFieldPolynomial sigma_i;
    while (true) {
        Zi = Zi_2 % Zi_1;
        GaloisFieldPolynomial qi = Zi_2 / Zi_1;
        sigma_i = sigma_i_2 - qi * sigma_i_1;

        // Check if the degree of remainder polynomial Zi is less the t = (n - k) / 2

        if (Zi.deg() < (n - k) / 2)
            break;

        // Update the polynomials accordingly
        Zi_2 = Zi_1;
        Zi_1 = Zi;

        sigma_i_2 = sigma_i_1;
        sigma_i_1 = sigma_i;

    }// while statement
    GaloisFieldElement sigma_i_Elements[n];
    for (int i = 0; i < n; i++) {
        sigma_i_Elements[i] = GaloisFieldElement(&galois_field, 0);
    }
    int numberOfError = 0;
    int *error_pos = new int[n];
    for (int i = 0; i < n; i++)
        error_pos[i] = 0;

    for (int i = 0; i <= n - 1; i++) {
        for (int j = 0; j <= sigma_i.deg(); j++) {
            sigma_i_Elements[n - i - 1] = sigma_i_Elements[n - i - 1] + galois_field.alpha(j * i) * sigma_i[j];
        }// for statement j
        if (sigma_i_Elements[n - i - 1] == 0) {
            error_pos[numberOfError] = n - i;
            numberOfError++; // find the number of errors detected
        }
    }// for statement i
    int errorPositions[numberOfError];
    for (int i = 0; i < numberOfError; ++i) {
        errorPositions[i] = error_pos[i];
    }
    delete[]error_pos; // deallocate memory

    GaloisFieldElement **beta = new GaloisFieldElement *[numberOfError];
    for (int i = 0; i < numberOfError; i++)
        beta[i] = new GaloisFieldElement[numberOfError];


    // Solve the linear system of equations to find the error values
    for (int i = 1; i <= numberOfError; i++) {
        for (int j = 0; j < numberOfError; j++) {
            beta[i - 1][j] = GaloisFieldElement(&galois_field,
                                                galois_field.alpha(i * (n - errorPositions[numberOfError - j - 1])));
        } // for statement j
    }// for statement i

    // Solve the linear equations using Cramer's rule 
    GaloisFieldElement errorValues[numberOfError];
    GaloisFieldElement denominatorDeterminant = determinant(beta, numberOfError);

    GaloisFieldElement **numerator = new GaloisFieldElement *[numberOfError];
    for (int i = 0; i < numberOfError; i++)
        numerator[i] = new GaloisFieldElement[numberOfError];

    for (int l = 0; l < numberOfError; l++) {

        for (int j = 0; j < numberOfError; j++) {
            for (int i = 0; i < numberOfError; i++) {
                numerator[i][j] = beta[i][j];
            }// for statement i
        }// for statement j

        for (int i = 0; i < numberOfError; i++) {
            numerator[i][l] = syndromeElements[n - k - i - 1];
        }// for statement i
        GaloisFieldElement numeratorDeterminant = determinant(numerator, numberOfError);
        errorValues[l] = numeratorDeterminant / denominatorDeterminant;
    }// for statement l

    GaloisFieldElement errors[n];
    for (int i = 0; i < n; i++)
        errors[i] = GaloisFieldElement(&galois_field, 0);

    for (int i = 0; i < numberOfError; i++)
        errors[n - errorPositions[numberOfError - i - 1]] = errorValues[i];

    // Correct for the error
    GaloisFieldPolynomial errorPolynomial = GaloisFieldPolynomial(&galois_field, n - 1, errors);
    cout << "encodedPolynomial: " << encodedPolynomial << endl;
    cout << "receivedPolynomial: " << receivedPolynomial << endl;

    cout << "errorPolynomial: " << errorPolynomial << endl;
    GaloisFieldPolynomial decodedPolynomial = receivedPolynomial + errorPolynomial;
    cout << "decodedPolynomial: " << decodedPolynomial << endl;
    GaloisFieldElement decodedMessage[k];
    for (int i = n - 1; i > n - k - 1; i--) {
        decodedMessage[i] = decodedPolynomial[i];
        cout << "decodedMessage: " << decodedMessage[i] << endl;
    }
    return 0;
}// end of main


GaloisFieldElement determinant(GaloisFieldElement **matrix, int n) {
    GaloisFieldElement det = GaloisFieldElement(&galois_field, 0);
    GaloisFieldElement **submatrix = new GaloisFieldElement *[n];
    for (int i = 0; i < n; i++)
        submatrix[i] = new GaloisFieldElement[n];


    if (n == 2)
        return ((matrix[0][0] * matrix[1][1]) + (matrix[1][0] * matrix[0][1]));
    else {
        for (int x = 0; x < n; x++) {
            int subi = 0;
            for (int i = 1; i < n; i++) {
                int subj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == x)
                        continue;
                    submatrix[subi][subj] = matrix[i][j];
                    subj++;
                }
                subi++;
            }

            det = det + (pow(1, x) * matrix[0][x] * determinant(submatrix, n - 1));
        }
    }
    return det;
}