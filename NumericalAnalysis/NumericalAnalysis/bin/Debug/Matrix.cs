namespace NumericSys
{
    using System;
    using System.Collections.Generic;

    //This structure will be reformed to eliminate the usage of any dynamic data type.
    public class Mat<T> : IEquatable<Mat<T>>, ICloneable where T : struct, IEquatable<T>, IComparable<T>
    {
        private T zero;
        private T[,] matrixArray;
        private int rows;
        private int columns;


        /// <summary>
        /// Default value of the struct data type.
        /// </summary>
        public T Zero => zero;
        public int Rows => rows;
        public int Columns => columns;
        /// <summary>
        /// Returning the matrix cell of i and j.
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <returns></returns>
        public T this[int i, int j]
        {
            get
            {
                return matrixArray[i, j];
            }
            set
            {
                matrixArray[i, j] = value;
            }
        }

        /// <summary>
        /// Factory matrix builder.
        /// </summary>
        public static MatrixFactory<T> Factory = new MatrixFactory<T>();


        /// <summary>
        /// 
        /// </summary>
        /// <param name="row"></param>
        /// <param name="column"></param>
        internal Mat(int rows, int columns)
        {
            zero = default(T);
            this.rows = rows;
            this.columns = columns;
            matrixArray = new T[rows, columns];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matrixArray[i, j] = zero;
                }
            }
        }


        internal Mat(T[,] array)
        {
            zero = default(T);
            rows = array.GetLength(0);
            columns = array.GetLength(1);
            var matArray = array.Clone() as T[,];
            this.matrixArray = matArray;

        }


        public bool Equals(Mat<T> other)
        {
            if (matrixArray.GetLength(0) != other.matrixArray.GetLength(0))
            {
                return false;
            }

            if (matrixArray.GetLength(1) != other.matrixArray.GetLength(1))
            {
                return false;
            }

            for (int i = 0; i < matrixArray.GetLength(0); i++)
            {
                for (int j = 0; j < matrixArray.GetLength(1); j++)
                {
                    if (!this[i, j].Equals(other.matrixArray[i, j]))
                    { return false; }

                }
            }

            return true;
        }


        public object Clone()
        {
            var mat = Factory.Create(rows, columns);

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    mat[i, j] = this[i, j];
                }
            }

            return mat;
        }


        public bool IsSquare()
        {
            return matrixArray.GetLength(0) == matrixArray.GetLength(1);
        }


        public List<T> Diagonals()
        {
            List<T> diag = new List<T>();
            for (int i = 0; i < rows - 1; i++)
            {
                diag.Add(this[i, i]);
            }
            return diag;
        }


        public virtual double ForbeniusNorm()
        {
            double totalSquareSum = 0;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    totalSquareSum += Multiply(matrixArray[i, j], matrixArray[i, j]);
                }
            }

            return Math.Sqrt(totalSquareSum);
        }


        public Mat<T> UpperTriangleMatrix()
        {
            var matrix = this.Clone() as Mat<T>;

            T temp = default(T);

            for (int i = 0; i < matrix.columns - 1; i++)
            {
                for (int j = (i + 1); j < matrix.columns; j++)
                {
                    for (int k = 0; k < matrix.columns; k++)
                    {
                        if (k == 0)
                        {
                            temp = matrix[j, i];
                        }
                        matrix[j, k] = SubtractType(matrix[j, k], (DivideType(MultiplyType(temp, matrix[i, k]), matrix[i, i])));
                    }
                }
                if (matrix[i + 1, i + 1].Equals(0))
                    break;
            }

            return matrix;
        }


        public double Determinant()
        {
            var mat = UpperTriangleMatrix();
            double accumliative = 1;
            for (int i = 0; i < rows; i++)
            {
                accumliative = MultiplyDouble(accumliative, mat[i, i]);
            }

            return accumliative;
        }


        public bool IsSingular()
        {
            return Determinant() == 0;
        }


        public Mat<T> Transpose()
        {
            var matResult = Factory.Create(columns, rows);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matResult[j, i] = matrixArray[i, j];
                }
            }
            return matResult;
        }


        public Mat<T> LowerTriangle()
        {
            return this.Multiply(UpperTriangleMatrix().Inverse());
        }


        public Mat<T> Inverse()
        {
            if (Determinant() == 0)
            {
                throw new InvalidOperationException("This matrix is a singular matrix");
            }
            if (!IsSquare())
            {
                throw new InvalidOperationException("This matrix is not a square matrix");
            }
            var matCofactor = Factory.Create(rows, columns);

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    dynamic x = Cofactor(i, j, rows).Determinant();
                    if ((i + j) % 2 == 0)
                        matCofactor[i, j] = x;
                    else
                        matCofactor[i, j] = -1 * x;
                }
            }

            var matTranspose = matCofactor.Transpose();
            dynamic scalar = Determinant();
            Mat<T> matResult = matTranspose.MutliplyScalar((1.00 / scalar));
            for (int i = 0; i < matResult.rows; i++)
            {
                for (int j = 0; j < matResult.columns; j++)
                {
                    dynamic x = matResult[i, j];
                    if (double.IsNaN(x))
                        matResult[i, j] = zero;
                }
            }
            return matResult;
        }


        public Mat<T> Multiply(Mat<T> matrix)
        {
            //Row of A.
            int rA = this.matrixArray.GetLength(0);
            //Column of A.
            int cA = this.matrixArray.GetLength(1);
            //Row of B.
            int rB = matrix.matrixArray.GetLength(0);
            //Column of B.
            int cB = matrix.matrixArray.GetLength(1);

            T temp = default(T);

            T[,] kHasil = new T[rA, cB];

            //Check if the matrices can be multiplied.
            if (cA != rB)
            {
                throw new Exception("Invalid Matrix Multiplication");
            }

            else
            {
                //The three loops to multiply the matrices.
                for (int i = 0; i < rA; i++)
                {
                    for (int j = 0; j < cB; j++)
                    {
                        temp = default(T);
                        for (int k = 0; k < cA; k++)
                        {
                            temp = AddType(temp, MultiplyType(this[i, k], matrix[k, j]));
                        }
                        kHasil[i, j] = temp;
                    }
                }
                return Factory.Create(kHasil);
            }
        }


        public Mat<T> MutliplyScalar(T scalar)
        {
            var matResult = this.Clone() as Mat<T>;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matResult[i, j] = MultiplyType(matrixArray[i, j], scalar);
                }
            }
            return matResult;
        }


        public void Negate()
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matrixArray[i, j] = Multiply(-1, matrixArray[i, j]);
                }
            }

        }


        public void Normalize()
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matrixArray[i, j] = Divide(matrixArray[i, j], ForbeniusNorm());
                }
            }
        }


        public Mat<T> Add(Mat<T> matrix)
        {
            var matResult = this.Clone() as Mat<T>;

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matResult[i, j] = AddType(matrixArray[i, j], matrix.matrixArray[i, j]);
                }
            }
            return matResult;
        }


        public Mat<T> Subtract(Mat<T> matrix)
        {
            var matResult = this.Clone() as Mat<T>;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    matResult[i, j] = SubtractType(matrixArray[i, j], matrix.matrixArray[i, j]);
                }
            }
            return matResult;
        }


        private static T Divide(T first, double second)
        {
            dynamic d1 = first;
            dynamic d2 = second;
            return (d1 / d2);
        }


        private static double MultiplyDouble(double scalar, T second)
        {
            dynamic d1 = scalar;
            dynamic d2 = second;
            return (d1 * d2);
        }


        private static double Add(T first, T second)
        {
            dynamic d1 = first;
            dynamic d2 = second;
            return (d1 + d2);
        }


        private static T AddType(T first, T second)
        {
            dynamic d1 = first;
            dynamic d2 = second;
            return (d1 + d2);
        }


        private static double Subtract(T first, T second)
        {
            dynamic d1 = first;
            dynamic d2 = second;
            return (d1 - d2);
        }


        private static T Subtract(double scalar, T second)
        {
            dynamic d1 = scalar;
            dynamic d2 = second;
            return (d1 - d2);
        }


        private static T Add(double scalar, T second)
        {
            dynamic d1 = scalar;
            dynamic d2 = second;
            return (d1 + d2);
        }


        private static double Multiply(T first, T second)
        {
            dynamic d1 = first;
            dynamic d2 = second;
            return (d1 * d2);
        }


        private static T MultiplyType(T first, T second)
        {
            dynamic d1 = first;
            dynamic d2 = second;
            return (d1 * d2);
        }


        private static T SubtractType(T first, T second)
        {
            dynamic d1 = first;
            dynamic d2 = second;
            return (d1 - d2);
        }


        private static T DivideType(T first, T second)
        {
            dynamic d1 = first;
            dynamic d2 = second;
            return (d1 / d2);
        }


        private static T Multiply(double scalar, T second)
        {
            dynamic d1 = scalar;
            dynamic d2 = second;
            return (d1 * d2);
        }


        private Mat<T> Cofactor(int row, int column, int size)
        {
            T[,] matArray = new T[size - 1, size - 1];

            List<T> matArray1D = new List<T>(size * size);

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    if (i == row || j == column)
                    {
                        continue;
                    }

                    matArray1D.Add(this[i, j]);
                }
            }
            int counter = 0;

            for (int i = 0; i < size - 1; i++)
            {
                for (int j = 0; j < size - 1; j++)
                {

                    matArray[i, j] = matArray1D[counter];

                    counter++;
                }

            }
            return Factory.Create(matArray);
        }

    }
}
