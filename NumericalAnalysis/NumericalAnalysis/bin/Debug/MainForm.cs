using System;
using System.Collections.Generic;
using System.Data;
using System.Windows.Forms;

namespace NumericSys
{
    public partial class MainForm : Form
    {
        DataTable tableLeftHand;
        DataTable tableRightHand;
        DataTable initialGuessTable;
        Mat<double> leftHandMatrix;
        Mat<double> rightHandMatrix;
        bool showFlag;
        int rows;
        int cols;

        public MainForm()
        {
            InitializeComponent();
            rows = 1;
            cols = 0;
            tableLeftHand = new DataTable();
            tableRightHand = new DataTable();
            initialGuessTable = new DataTable();
            dataGridView1.DataSource = tableLeftHand;
            dataGridView2.DataSource = tableRightHand;
            dataGridViewInitGuess.DataSource = initialGuessTable;
            showFlag = true;
        }

        private void button1_Click_1(object sender, EventArgs e)
        {
            tableLeftHand.Rows.Add();
            rows++;
        }

        private void button2_Click(object sender, EventArgs e)
        {
            tableLeftHand.Columns.Add($"{cols + 1}");
            cols++;
        }

        private void btnInitSol_Click(object sender, EventArgs e)
        {
            tableRightHand.Columns.Add();
            for (int i = 0; i < rows-1; i++)
            {
                tableRightHand.Rows.Add();
            }
        }

        private void FillMatrices()
        {
            leftHandMatrix = Mat<double>.Factory.Create(rows, cols);
            rightHandMatrix = Mat<double>.Factory.Create(rows, 1);
            for (int i = 0; i < tableLeftHand.Rows.Count; i++)
                for (int j = 0; j < tableLeftHand.Columns.Count; j++)
                {
                    object o = tableLeftHand.Rows[i].ItemArray[j];
                    string oString = o.ToString();
                    double value = int.Parse(oString);
                    leftHandMatrix[i, j] = value;
                }

            for (int i = 0; i < tableRightHand.Rows.Count; i++)
            {
                object o = tableRightHand.Rows[i].ItemArray[0];
                string oString = o.ToString();
                double value = int.Parse(oString);
                rightHandMatrix[i, 0] = value;
            }
        }

        private void button4_Click(object sender, EventArgs e)
        {
            FillMatrices();

            List<double> doubleList = new List<double>();
            for (int i = 0; i < leftHandMatrix.Columns; i++)
            {
                doubleList.Add(leftHandMatrix[0, i]);
            }

            double root = NumericalUtils.Bisection(double.Parse(txtBoxInitBis.Text), double.Parse(txtBoxSecBis.Text), double.Parse(txtBoxTol.Text), int.Parse(txtBoxIter.Text), doubleList, NumericalUtils.CalcF, out int counter);

            txtBoxRootBis.Text = root.ToString();
            txtBoxNitBis.Text = counter.ToString();
        }

        private void button5_Click(object sender, EventArgs e)
        {
            FillMatrices();

            List<double> doubleList = new List<double>();
            for (int i = 0; i < leftHandMatrix.Columns; i++)
            {
                doubleList.Add(leftHandMatrix[0, i]);
            }

            double root = NumericalUtils.FalsePosition(double.Parse(txtBoxInitFls.Text), double.Parse(txtBoxSecFls.Text), double.Parse(txtBoxTol.Text), int.Parse(txtBoxIter.Text), doubleList, NumericalUtils.CalcF, out int counter);

            txtBoxRootFls.Text = root.ToString();
            txtBoxNitFls.Text = counter.ToString();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            FillMatrices();

            List<double> doubleList = new List<double>();
            for (int i = 0; i < leftHandMatrix.Columns; i++)
            {
                doubleList.Add(leftHandMatrix[0, i]);
            }

            var area = NumericalUtils.TrapezoidalMethod(double.Parse(txtBoxLowerTrap.Text), double.Parse(txtBoxUpperTrap.Text), int.Parse(txtBoxNoStripsTrap.Text), doubleList, NumericalUtils.CalcF);

            txtBoxAreaTrap.Text = area.ToString();

        }

        private void button6_Click(object sender, EventArgs e)
        {
            FillMatrices();

            List<double> doubleList = new List<double>();
            for (int i = 0; i < leftHandMatrix.Columns; i++)
            {
                doubleList.Add(leftHandMatrix[0, i]);
            }

            var area = NumericalUtils.MidPointMethod(double.Parse(txtBoxLowerMid.Text), double.Parse(txtBoxUpperMid.Text), int.Parse(txtBoxNoStripsMid.Text), doubleList, NumericalUtils.CalcF);

            txtBoxAreaMid.Text = area.ToString();
        }

        private void button7_Click(object sender, EventArgs e)
        {
            FillMatrices();

            var maxEigen = NumericalUtils.GetMaxEigen(leftHandMatrix, out double[] eigenVector, int.Parse(txtBoxIter.Text), double.Parse(txtBoxTol.Text), out int counter);
            DataTable table = new DataTable();
            table.Columns.Add("1");
            for (int i = 0; i < eigenVector.Length; i++)
            {
                table.Rows.Add(eigenVector[i]);
            }
            dataGridViewEigen.DataSource = table;
            dataGridViewEigen.Enabled = false;
            txtBoxMaxEigen.Text = maxEigen.ToString();
            txtBoxiterPower.Text = counter.ToString();
        }

        private void groupBox6_Enter(object sender, EventArgs e)
        {

        }

        private void button8_Click(object sender, EventArgs e)
        {
            FillMatrices();
            Mat<double> initMat = Mat<double>.Factory.Create(rightHandMatrix.Rows, rightHandMatrix.Columns);
            for (int i = 0; i < initialGuessTable.Rows.Count; i++)
            {
                object o = initialGuessTable.Rows[i].ItemArray[0];
                string oString = o.ToString();
                double value = int.Parse(oString);
                initMat[i, 0] = value;
            }
            var approxMat = NumericalUtils.GaussIterative(leftHandMatrix, rightHandMatrix, initMat, double.Parse(txtBoxTol.Text), int.Parse(txtBoxIter.Text), out int counter);
            double[] arr = new double[approxMat.Rows];

            for (int i = 0; i < approxMat.Rows; i++)
            {
                arr[i] = approxMat[i, 0];
            }

            DataTable table = new DataTable();
            table.Columns.Add("1");
            for (int i = 0; i < arr.Length; i++)
            {
                table.Rows.Add(arr[i]);
            }
            dataGridViewApproxRes.DataSource = table;
            dataGridViewApproxRes.Enabled = false;

            txtBoxGaussSidel.Text = counter.ToString();


        }

        private void button9_Click(object sender, EventArgs e)
        {
            if (showFlag)
            {
                FillMatrices();

                initialGuessTable.Columns.Add("1");

                for (int i = 0; i < rightHandMatrix.Rows - 1; i++)
                {
                    initialGuessTable.Rows.Add();
                }
                showFlag = false;
            }
        }

        private void button10_Click(object sender, EventArgs e)
        {
            FillMatrices();

            var maxEigen = NumericalUtils.InversePowerIteration(leftHandMatrix, out double[] eigenVector, int.Parse(txtBoxIter.Text), double.Parse(txtBoxTol.Text), out int counter);
            DataTable table = new DataTable();
            table.Columns.Add("1");
            for (int i = 0; i < eigenVector.Length; i++)
            {
                table.Rows.Add(eigenVector[i]);
            }
            dataGridViewInv.DataSource = table;
            dataGridViewInv.Enabled = false;
            txtBoxMinEigen.Text = maxEigen.ToString();
            txtBoxiterInvPower.Text = counter.ToString();
        }
    }
}
