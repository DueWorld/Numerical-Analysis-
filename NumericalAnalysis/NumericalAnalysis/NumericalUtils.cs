namespace NumericSys
{
    using System;
    using System.Collections.Generic;

    public static class NumericalUtils
    {

       /// <summary>
       /// Calculating a polynomial function, serves as a delegate for any polynomial serving algorithm.
       /// </summary>
       /// <param name="xValue"></param>
       /// <param name="factors"></param>
       /// <returns></returns>
        public static double CalcF(double xValue, List<double> factors)
        {
            double sum = 0;

            double power = factors.Count - 1;

            for (int i = 0; i < factors.Count; i++)
            {
                var test = Math.Pow(xValue, power) * factors[i];
                sum += Math.Pow(xValue, power) * factors[i];
                power--;
            }

            return sum;
        }




        /// <summary>
        /// Gets the integration (Area under curve) using trapezoidal method technique.
        /// </summary>
        /// <param name="lower"></param>
        /// <param name="upper"></param>
        /// <param name="noStrip"></param>
        /// <param name="factors"></param>
        /// <param name="polyFunction"></param>
        /// <returns></returns>
        public static double TrapezoidalMethod(double lower, double upper, int noStrip, List<double> factors, Func<double, List<double>, double> polyFunction)
        {
            //starting initial of sum.
            double sum = 0;

            List<double> vecX = new List<double>(noStrip);

            List<double> vecY = new List<double>(noStrip);

            //Calculating the delta X.
            double deltaX = (upper - lower) / (noStrip - 1);

            for (int i = 0; i < noStrip; i++)
            {
                vecX.Add(lower + (i) * deltaX);
            }

            foreach (double x in vecX)
            {
                vecY.Add(polyFunction(x, factors));
            }

            //Tracking the repetitions of the f(x).
            for (int i = 0; i < noStrip; i++)
            {
                if(i==0 || i== noStrip-1)
                {
                    sum += vecY[i] / 2;
                }
                else
                {
                    sum += vecY[i];
                }
            }

            //After common factor.

            return sum * deltaX;
        }


        /// <summary>
        /// Gets the integration (Area under curve) using mid point method technique.
        /// </summary>
        /// <param name="lower"></param>
        /// <param name="upper"></param>
        /// <param name="noStrip"></param>
        /// <param name="factors"></param>
        /// <param name="polyFunction"></param>
        /// <returns></returns>
        public static double MidPointMethod(double lower, double upper, int noStrip, List<double> factors, Func<double, List<double>, double> polyFunction)
        {
            //Calculating the delta x value.
            double deltaX = (upper - lower) / noStrip;

            //Getting the mid x values.
            double xValue = lower + deltaX * 0.5;
            //Initiating the sum.
            double sum = 0;

            for (int i = 0; i < noStrip; i++)
            {
                sum += polyFunction(xValue, factors);
                xValue += deltaX;
            }
            //Taking back the common factor of the delta x.
            return sum * deltaX;
        }



        /// <summary>
        /// Performs a gauss iterative on a linear system to get the approximate solution.
        /// </summary>
        /// <param name="leftHandSide"></param>
        /// <param name="rightHandSide"></param>
        /// <param name="initGuess"></param>
        /// <param name="tolerance"></param>
        /// <param name="iterationNumber"></param>
        /// <param name="counter"></param>
        /// <returns></returns>
        public static Mat<double> GaussIterative(Mat<double> leftHandSide, Mat<double> rightHandSide,Mat<double> initGuess, double tolerance, int iterationNumber, out int counter)
        {
            double eTolCheck = tolerance + 1;

            counter = 0;

            double auxilary = 0;

            double segma = 0;

            while( eTolCheck> tolerance && counter<iterationNumber)
            {
                if(counter > 0)
                {
                    auxilary = initGuess[0,0];
                }

                for (int i = 0; i < leftHandSide.Rows; i++)
                {
                    segma = 0;
                    for (int j = 0; j < leftHandSide.Columns; j++)
                    {
                        if(j!= i)
                        {
                            segma += leftHandSide[i, j] * initGuess[j, 0];
                        }
                    }
                        initGuess[i, 0] = (1 / leftHandSide[i, i]) * (rightHandSide[i, 0] - segma);
                }

                if(counter==0)
                {
                    auxilary = initGuess[0, 0];
                }
                else
                {
                    eTolCheck = Math.Abs((initGuess[0, 0] - auxilary) / initGuess[0, 0]);
                }

                counter++;


            }

            return initGuess;
        }


        /// <summary>
        /// Getting the maximum Eigen values with the power iteration method.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="eigenVector"></param>
        /// <returns></returns>
        public static double GetMaxEigen(Mat<double> A, out double[] eigenVector, int iterationNumber, double epsilon, out int counter)
        {
            counter = 1;
            //Row size of the matrix.
            int rowSize = A.Rows;
            //Row size of the matrix.
            int columnSize = A.Columns;

            //The scaled image vector X.
            Mat<double> x = Mat<double>.Factory.Create(rowSize, 1);
            //Checks if this matrix is a square matrix or not.
            if (rowSize != columnSize)
            {
                throw new Exception("This matrix is not a square matrix");
            }
            //Filling the x with an initial value of 1
            //so the X will be [1,1,1...].
            for (int i = 0; i < rowSize; i++)
            {
                x[i, 0] = 1;
            }

            //Initializing the new and old norm variables.
            double newNorm = 0;

            double oldNorm = 0;

            //Imaging the X (This algorithm didn't use another variable (B) in slided.
            //Just the X being modified and reiterated upon.
            x = A.Multiply(x);

            //Calculating the norm.
            newNorm = x.ForbeniusNorm();

            //Scaling the vector X with 1/lambda.
            x = x.MutliplyScalar(1 / newNorm);

            //The main loop with the K(max) iterations limit.
            for (int i = 0; i < iterationNumber; i++)
            {
                //Equating the old norm by the new norm
                //As the new norm will be changed in the next few lines.
                oldNorm = newNorm;

                //Re imaging the vector.
                x = A.Multiply(x);

                //The new norm updated for the re-imaged vector.
                newNorm = x.ForbeniusNorm();

                //ReScaling the x.
                x = x.MutliplyScalar(1 / newNorm);

                //Convergence condition, also taking care of the absolute value.
                if (Math.Abs(((newNorm - oldNorm) / newNorm)) <= epsilon)
                {
                    break;
                }
                //If the loop continues the X will be re-imaged then rescaled,
                //and new values for the newNorm and oldNorm will appear.
                counter++;
            }

            //Populating the Eigen vector.
            double[] temp = new double[rowSize];

            for (int i = 0; i < rowSize; i++)
            {
                temp[i] = x[i, 0];
            }

            //Returning the Eigen vector.
            eigenVector = temp;

            //Returning the last norm which will be the maximum Eigen value.
            return newNorm;
        }

        /// <summary>
        /// Perform an inverse power iteration method on a square matrix.
        /// </summary>
        /// <param name="A"></param>
        /// <param name="eigenVector"></param>
        /// <param name="iterationNumber"></param>
        /// <param name="epsilon"></param>
        /// <param name="counter"></param>
        /// <returns></returns>
        public static double InversePowerIteration(Mat<double> A, out double[] eigenVector, int iterationNumber, double epsilon, out int counter)
        {
            counter = 1;
            //Row size of the matrix.
            int rowSize = A.Rows;
            //Row size of the matrix.
            int columnSize = A.Columns;

            //The scaled image vector X.
            Mat<double> x = Mat<double>.Factory.Create(rowSize, 1);
            //Checks if this matrix is a square matrix or not.
            if (rowSize != columnSize)
            {
                throw new Exception("This matrix is not a square matrix");
            }
            //Filling the x with an initial value of 1
            //so the X will be [1,1,1...].
            for (int i = 0; i < rowSize; i++)
            {
                x[i, 0] = 1;
            }

            //Initializing the new and old norm variables.
            double newNorm = 0;

            double oldNorm = 0;

            var matInverse = A.Inverse();

            //Imaging the X (This algorithm didn't use another variable (B) in slided.
            //Just the X being modified and reiterated upon.
            x = matInverse.Multiply(x);

            //Calculating the norm.
            newNorm = x.ForbeniusNorm();

            //Scaling the vector X with 1/lambda.
            x = x.MutliplyScalar(1 / newNorm);

            //The main loop with the K(max) iterations limit.
            for (int i = 0; i < iterationNumber; i++)
            {
                //Equating the old norm by the new norm
                //As the new norm will be changed in the next few lines.
                oldNorm = newNorm;

                //Re imaging the vector.
                x = matInverse.Multiply(x);

                //The new norm updated for the re-imaged vector.
                newNorm = x.ForbeniusNorm();

                //ReScaling the x.
                x = x.MutliplyScalar(1 / newNorm);

                //Convergence condition, also taking care of the absolute value.
                if (Math.Abs(((newNorm - oldNorm) / newNorm)) <= epsilon)
                {
                    break;
                }
                //If the loop continues the X will be re-imaged then rescaled,
                //and new values for the newNorm and oldNorm will appear.
                counter++;
            }

            //Populating the Eigen vector.
            double[] temp = new double[rowSize];

            for (int i = 0; i < rowSize; i++)
            {
                temp[i] = x[i, 0];
            }

            //Returning the Eigen vector.
            eigenVector = temp;

            //Returning the last norm which will be the maximum Eigen value.
            return (1/newNorm);

        }


        /// <summary>
        /// Using bisection open method for a polynomial function to find its root.
        /// </summary>
        /// <param name="xInitial"></param>
        /// <param name="xSecondary"></param>
        /// <param name="tolerance"></param>
        /// <param name="iterationNo"></param>
        /// <param name="factors"></param>
        /// <param name="polyFunction"></param>
        /// <param name="counter"></param>
        /// <returns></returns>
        public static double Bisection(double xInitial, double xSecondary, double tolerance, int iterationNo, List<double> factors, Func<double, List<double>, double> polyFunction, out int counter)
        {
            //Initiating the initial and secondary functions and x values.
            double yInitial = polyFunction(xInitial, factors);
            counter = 1;
            double ySecondary = polyFunction(xSecondary, factors);

            //Check if any is already the root.
            if (yInitial == 0)
                return xInitial;
            if (ySecondary == 0)
                return xSecondary;
            if (yInitial * ySecondary > 0)
                throw new Exception("No roots are found");

            //performs the first run of the initial and secondary bisecting.
            double xMean = (xInitial + xSecondary) / 2;
            double yMean = polyFunction(xMean, factors);
            double eTolCheck = 1 + tolerance;
            double auxilaryVariable = 0;

            //The main while loop with the tolerance and counter checks.
            while (yMean != 0 && eTolCheck > tolerance && counter < iterationNo)
            {
                //Not to get two intervals both having the same sign, so it goes away from the root.
                if (yInitial * yMean < 0)
                {
                    xSecondary = xMean;
                    ySecondary = yMean;
                }
                else
                {
                    xInitial = xMean;
                    yInitial = yMean;
                }

                //Update the mean values.
                auxilaryVariable = xMean;
                xMean = (xInitial + xSecondary) / 2;
                yMean = polyFunction(xMean, factors);

                //Update the current tolerance.
                eTolCheck = Math.Abs((xMean - auxilaryVariable) / xMean);
                counter++;
            }
            return xMean;
        }


        /// <summary>
        /// Using false position open method for a polynomial function to find its root.
        /// </summary>
        /// <param name="xInitial"></param>
        /// <param name="xSecondary"></param>
        /// <param name="tolerance"></param>
        /// <param name="iterationNo"></param>
        /// <param name="factors"></param>
        /// <param name="polyFunction"></param>
        /// <param name="counter"></param>
        /// <returns></returns>
        public static double FalsePosition(double xInitial, double xSecondary, double tolerance, int iterationNo, List<double> factors, Func<double, List<double>, double> polyFunction, out int counter)
        {
            //Initiating the initial and secondary functions and x values.
            double yInitial = polyFunction(xInitial, factors);

            counter = 1;
            double eTolCheck = 1 + tolerance;
            double auxilaryVariable = 0;
            double ySecondary = polyFunction(xSecondary, factors);

            //Check if any is already the root.

            if (yInitial == 0)
                return xInitial;
            if (ySecondary == 0)
                return xSecondary;
            if (yInitial * ySecondary > 0)
                throw new Exception("No roots are found");
            double slope = (yInitial - ySecondary) / (xInitial - xSecondary);

            //performs the first run of the initial and secondary bisecting.
            double xMean = ((-yInitial) + slope * xInitial) / slope;
            double yMean = polyFunction(xMean, factors);

            //The main while loop with the tolerance and counter checks.
            while (yMean != 0 && eTolCheck > tolerance && counter < iterationNo)
            {
                //Not to get two intervals both having the same sign, so it goes away from the root.

                if (yInitial * yMean < 0)
                {
                    xSecondary = xMean;
                    ySecondary = yMean;
                }
                else
                {
                    xInitial = xMean;
                    yInitial = yMean;
                }

                //Update the mean values.
                auxilaryVariable = xMean;
                slope = (yInitial - ySecondary) / (xInitial - xSecondary);
                xMean = ((-yInitial) + slope * xInitial) / slope;
                yMean = polyFunction(xMean, factors);

                //Update the current tolerance.
                eTolCheck = Math.Abs((xMean - auxilaryVariable) / xMean);
                counter++;
            }
            return xMean;
        }

    }
}
