namespace NumericSys
{
    using System;

    public class MatrixFactory<T> where T : struct, IEquatable<T>, IComparable<T>
    {

        internal MatrixFactory()
        {

        }

        public Mat<T> Create(int rows, int cols)
        {
            return new Mat<T>(rows, cols);
        }

        public Mat<T> Create(T[,] array)
        {
            return new Mat<T>(array);
        }
    }
}
