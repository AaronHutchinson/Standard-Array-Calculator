#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <windows.h>

using namespace std;

int vectorDimension = 0;
int codeDimension = 0;
int fieldChar = 0;



std::vector<int> convert( int integerToConvert )                            // Converts the given integer into base fieldChar and returns a vector of length codeDimension
{
    static std::vector<int> coefVecToReturn(codeDimension);
    for ( int i = 0 ; i < codeDimension ; i++ )
    {
        coefVecToReturn[i] = integerToConvert / pow( fieldChar, codeDimension - i - 1);
        integerToConvert = integerToConvert - coefVecToReturn[i] * pow(fieldChar, codeDimension - i - 1);

    }
    return coefVecToReturn;
}


std::vector<int> convertVectorDimension( int integerToConvert )                            // Converts the given integer into base fieldChar and returns a vector of length vectorDimension
{

    static std::vector<int> coefVecToReturn(vectorDimension);
    for ( int i = 0 ; i < vectorDimension ; i++ )
    {
        coefVecToReturn[i] = integerToConvert / pow( fieldChar, vectorDimension - i - 1);
        integerToConvert = integerToConvert - coefVecToReturn[i] * pow(fieldChar, vectorDimension - i - 1);

    }
    return coefVecToReturn;
}


std::vector<int> vectorAddition( std::vector<int> vec1 , std::vector<int> vec2 )
{
    std::vector<int> vectorSum(vectorDimension);
    for (int i = 0 ; i < vectorDimension ; i++ )
    {
        vectorSum[i] = ( vec1[i] + vec2[i] ) % fieldChar;
    }
    return vectorSum;
}

std::vector<int> linearCombination( std::vector<int> coefficientVector , std::vector< std::vector<int> > basis )            // Computes the linear combination of the basis vectors using the coefficient vector
{
    std::vector<int> combinationVector(vectorDimension,0);
    for ( int i = 0 ; i < vectorDimension ; i++ )
    {
        for ( int j = 0 ; j < basis.size() ; j++ )
        {
            combinationVector[i] = combinationVector[i] + coefficientVector[j] * basis[j][i];
        }
        combinationVector[i] = combinationVector[i] % fieldChar;
    }
    return combinationVector;
}

int weight( std::vector<int> vec)                          // Computes the weight of a vector
{
    int weightAccum = 0;
    for ( int i = 0 ; i < vectorDimension ; i++ )
    {
        if (vec[i] != 0)
            weightAccum++;
    }
    return weightAccum;
}

int vectorComparison( std::vector<int> vec1, std::vector<int> vec2)             // Compares two vectors of the same length to see if they are equal. Return 1 if they are equal, 0 if not
{
    int parityCheck = 1;
    for( int i = 0 ; i < vectorDimension ; i++ )
    {
        if ( vec1[i] != vec2[i] )
        {
            parityCheck = 0;
            break;
        }
    }
    return parityCheck;
}


std::vector< std::vector< std::vector<int> > > computeArray( std::vector< std::vector<int> > code, std::vector< std::vector<int> > vectorSpace )
{

    std::vector< vector< vector<int> > > standardArray( (int) pow( fieldChar, vectorDimension - codeDimension) );      // This vector holds the entire array. Think of it as a matrix: the first index is the row, second index is the column, third is the coordinate of the vector in row i and column j
    for( int i = 0 ; i < (int) pow( fieldChar, vectorDimension - codeDimension ) ; i++ )                             // Resize the array
    {
        standardArray[i].resize( (int) pow( fieldChar, codeDimension ) );
        for( int j = 0 ; j < (int) pow( fieldChar, codeDimension ) ; j++ )
        {
            standardArray[i][j].resize( vectorDimension );
            if ( i == 0 )                                                                                       // Fill in the first row of the array as the code
            {
                standardArray[i][j] = code[j];
            }
        }
    }


    int rowTracker = 1;                                                                                         // Keeps track of how many cosets have already been computed.
    int parityCheck = 0;                                                                                        // 1 if the vector we're looking for shows up in the array already, 0 if not
    std::vector<int> cosetRep(vectorDimension);                                                                 // Stores a coset representative
    while( rowTracker < (int) pow( fieldChar, vectorDimension - codeDimension ))
    {
        for( int i = 0 ; i < (int) pow( fieldChar, vectorDimension) ; i++ )
        {
            for( int j = 0 ; j < rowTracker ; j++ )
            {
                for( int k = 0 ; k < (int) pow( fieldChar, codeDimension) ; k++ )
                {
                    if ( vectorComparison( vectorSpace[i] , standardArray[j][k] ) == 1 )
                    {
                        parityCheck = 1;
                        break;
                    }
                }
                if ( parityCheck == 1)
                    break;
            }
            if ( parityCheck == 0 )                                                                             // If we found a vector which is not in the array already, add its coset to the array
            {
                cosetRep = vectorSpace[i];
                for( int j = 1 ; j < (int) pow( fieldChar, codeDimension ) ; j++)                                     // Find the minimal weight vector in this coset and store it in cosetRep
                {
                    if ( weight(vectorAddition( standardArray[0][j] , vectorSpace[i] )) < weight( cosetRep ))
                        cosetRep = vectorAddition( standardArray[0][j] , vectorSpace[i] );
                }
                for( int j = 0 ; j < (int) pow( fieldChar, codeDimension) ; j++ )                                     // Add the coset to the array with the minimal weight vector (cosetRep) as the first entry.
                {
                    standardArray[rowTracker][j] = vectorAddition( cosetRep, standardArray[0][j] );
                }
                rowTracker++;
            }
            parityCheck = 0;
        }
    }
    return standardArray;
}

std::vector< std::vector<int> > computeSpace()
{
    std::vector< std::vector<int> > standardBasis(vectorDimension);         // Define the standard basis
    for( int i = 0 ; i < vectorDimension ; i++ )
    {
        standardBasis[i].resize( vectorDimension );
        for ( int j = 0 ; j < vectorDimension ; j++)
        {
            if ( i == j)
                standardBasis[i][j] = 1;
            else
                standardBasis[i][j] = 0;
        }
    }


    // Now compute the entire vector space by taking the span of this basis
    std::vector< std::vector<int> > vectorSpace( pow( fieldChar, vectorDimension) );

    for( int i = 0 ; i < (int) pow(fieldChar , vectorDimension) ; i++ )
    {
        vectorSpace[i].resize(vectorDimension);
        vectorSpace[i] = linearCombination( convertVectorDimension(i), standardBasis );
    }
    return vectorSpace;

}

void printArray( std::vector< std::vector< std::vector<int> > > standardArray )
{
    for( int i = 0 ; i < pow( fieldChar, vectorDimension - codeDimension) ; i++ )
    {
        for( int j = 0 ; j < pow( fieldChar, codeDimension ) ; j++ )
        {
            for( int k = 0 ; k < vectorDimension ; k++ )
                cout << standardArray[i][j][k];
            cout << "  ";
        }
        cout << "\n   ";
    }
}

void printArrayToFile( std::vector< std::vector< std::vector<int> > > standardArray, std::vector< std::vector<int> > basis )
{
    ofstream outputFile;
    outputFile.open("Standard Array.txt");
    outputFile << "The following was generated using the Standard Array Calculator by Aaron Hutchinson with the data:\n   Vector Space Dimension: " << vectorDimension << "\n           Code Dimension: " << codeDimension << "\n     Field Characteristic: " << fieldChar << "\n   " << "                 Basis: ";

    for( int i = 0 ; i < codeDimension ; i++ )
    {
        for( int j = 0 ; j < vectorDimension ; j++ )
        {
            outputFile << basis[i][j];
        }
        outputFile << "\n                           ";
    }

    outputFile << "\nThe standard array for this data is given below. The first row of the array is the code, generated by the given basis. The remaining columns are the cosets of the code in the vector space, with the minimal weight vector of each coset displayed in the leftmost column.\n\n   ";
    for( int i = 0 ; i < pow( fieldChar, vectorDimension - codeDimension) ; i++ )
    {
        for( int j = 0 ; j < pow( fieldChar, codeDimension ) ; j++ )
        {
            for( int k = 0 ; k < vectorDimension ; k++ )
                outputFile << standardArray[i][j][k];
            outputFile << "  ";
        }
        outputFile << "\n   ";
    }
}

int main()
{
    // This program generates the standard array for a given linear code over Z_p.
    // The vector space dimension, code dimension, field characteristic, and code basis must be given.
    // With this information, the standard array is computed.



    // ***************************** Initialization Stage ************************
    cout << "We will compute the standard array for a given linear code. The field must be Z_p for some prime p.\nPlease enter the following information:\n";
    cout << "   Vector space dimension: n = ";
    cin >> vectorDimension;
    cout << "   Code dimension:         k = ";
    cin >> codeDimension;
    cout << "   Field characteristic:   p = ";
    cin >> fieldChar;
    cin.ignore();

    std::vector< std::vector<int> > basis(codeDimension);        // This vector will store the k-dimensional basis. The second dimension is the index in the basis, while the first dimension is the index of a digit in a given basis vector.
    for ( int i = 0 ; i < codeDimension ; i++ )                  // Initial the vector 'basis' to the correct size with all 0's
    {
        basis[i].resize(vectorDimension);
        for ( int j = 0 ; j < vectorDimension ; j++ )
        {
            basis[i][j] = 0;
        }
    }

    cout << "\nEnter " << codeDimension << " strings of length " << vectorDimension << " consisting of the integers 0 - " << fieldChar - 1 << " as the basis for the code. These vectors are assumed to be a basis; if they are not, there will be problems.\n   1. ";

    for (int i = 0; i < codeDimension ; i = i + 1)
    {
        string receivedInput;
        getline( cin, receivedInput );

        for ( int j = 0 ; j < vectorDimension ; j++ )
        {
            basis[i][j] = receivedInput[j] - '0';
        }

        if (i != codeDimension - 1)
        {
            cout << "   " << i+2 << ": ";
        }

    }


    // **************************** Compute the code *****************************
    // We will now use the basis to compute the code in the array 'code'.
    double codeSize = pow(fieldChar,codeDimension);
    std::vector<vector<int> > code(codeSize);               // This 'array' will hold all vectors in the code
    for ( int i = 0 ; i < codeSize ; i++ )                  // Initial the vector 'code' to the correct size with all 0's
        code[i].resize(vectorDimension);

    std::vector<int> coefficientVector(codeDimension);      // This array will hold the coefficients of a vector defined in terms of the basis vectors
    for ( int i = 0 ; i < codeDimension ; i++ )             // Initial the coefficient vector to 0
            coefficientVector[i] = 0;


    cout << "\n\n";
    for ( int i = 0 ; i < codeSize ; i++ )                  // Compute all vectors in the code
    {
        coefficientVector = convert(i);
        code[i] = linearCombination( coefficientVector, basis );
    }


    std::vector< std::vector<int> > vectorSpace( pow( fieldChar, vectorDimension ) );
    vectorSpace = computeSpace();
    std::vector< std::vector< std::vector<int> > > standardArray = computeArray( code, vectorSpace );

    cout << "The standard array is given below. The first row of the array is the code, generated by the given basis. The remaining columns are the cosets of the code in the vector space, with the minimal weight vector of each coset displayed in the leftmost column. A .txt file has also been created which contains this data.\n\n   ";

    if ( 3 + vectorDimension * pow( fieldChar, codeDimension )  + 2*(codeSize - 1) < 100 ) // If the array is small enough to fit in the command prompt, the program will display it there. If not, a .txt file will be generated.
    {
        printArray( standardArray );
        cout << "\n\n\n";
    }
    else
        cout << "The array is too big to display in command prompt. Please refer to the .txt file.\n\n\n";


    printArrayToFile( standardArray, basis );

    return 0;
}
