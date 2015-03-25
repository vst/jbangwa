/*
 * Copyright (c) 2015 Vehbi Sinan Tunalioglu <vst@vsthost.com>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.vsthost.rnd.jbangwa;

import com.vsthost.rnd.commons.math.ext.linear.DMatrixUtils;
import com.vsthost.rnd.commons.math.ext.linear.EMatrixUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * Provides convenience functionalities for finance-related computational tasks.
 */
public class FinanceUtils {

    /**
     * Computes the simple returns of the prices provided.
     *
     * @param prices The prices of which the returns are going to be calculated.
     * @return The return series as a double array.
     */
    public static double[] getReturns (final double[] prices) {
        // Declare and initialize the return value:
        double[] returns = new double[prices.length - 1];

        // Iterate over prices and calculate returns:
        // (Indexing is done by "i + 1" instead of more intuitive "i - 1"
        // to avoid redundancy).
        for (int i = 0; i < returns.length; i++) {
            returns[i] = (prices[i + 1] - prices[i]) / prices[i];
        }

        // Done, return returns:
        return returns;
    }

    /**
     * Computes the simple returns of the prices provided.
     *
     * @param prices A list of prices where each column refers to a price series of a single asset.
     * @return The list of returns where each column refers to a return series of a single asset.
     */
    public static double[][] getReturns (final double[][] prices) {
        // Check prices:
        if (prices == null || prices.length == 0) {
            return prices;
        }

        // Compute dimensions:
        int rowCount = prices.length - 1;
        int colCount = prices[0].length;

        // Declare the return value:
        double[][] returns = new double[rowCount][colCount];

        // Iterate over price series and compute each return series:
        for (int row = 0; row < rowCount; row++) {
            for (int col = 0; col < colCount; col++) {
                returns[row][col] = (prices[row + 1][col] - prices[row][col]) / prices[row][col];
            }
        }

        // Done, return returns:
        return returns;
    }

    /**
     * Convenience function for return series calculation for RealMatrix based price series.
     *
     * @param prices A matrix of prices.
     * @return A matrix of return series.
     */
    public static RealMatrix getReturns (final RealMatrix prices) {
        return org.apache.commons.math3.linear.MatrixUtils.createRealMatrix(FinanceUtils.getReturns(prices.getData()));
    }

    /**
     * Detrends return series.
     *
     * @param returns Return series to be detrended.
     * @return Detrended return series.
     */
    public static RealMatrix detrendReturns(final RealMatrix returns) {
        // Subtract col means from rows and return:
        return EMatrixUtils.rowSubtract(returns, EMatrixUtils.colMeans(returns));
    }

    /**
     * Computes the returns of a portfolio with given returns and weights per asset.
     *
     * @param returns Per asset returns.
     * @param weights Per asset wright.
     * @return The portfolio return.
     */
    public static double computeReturns (double[] returns, double[] weights) {
        // Initialize the return value:
        double total = 0.0;

        // Iterate over the assets and add to total:
        for (int i = 0; i < returns.length; i++) {
            total += returns[i] * weights[i];
        }

        // Done, return the total return:
        return total;
    }

    /**
     * Computes portfolio returns of a series of asset returns and weights vector.
     *
     * @param returns Asset returns.
     * @param weights Asset weights.
     * @return Portfolio returns.
     */
    public static double[] portfolioReturns(final RealMatrix returns, final double[] weights) {
        // Declare the return value:
        double[] retval = new double[returns.getRowDimension()];

        // Iterate over rows:
        for (int row = 0; row < returns.getRowDimension(); row++) {
            for (int col = 0; col < returns.getColumnDimension(); col++) {
                retval[row] += returns.getEntry(row, col) * weights[col];
            }
        }

        // Done, return portfolio returns:
        return retval;
    }

    /**
     * Computes equal weighted return series.
     *
     * @param returns Asset returns.
     * @return Portfolio returns.
     */
    public static double[] portfolioReturns(final RealMatrix returns) {
        // Compute equal weights and call overload:
        return FinanceUtils.portfolioReturns(returns, FinanceUtils.equalWeights(returns.getColumnDimension()));
    }

    /**
     * Computes a vector of equal asset weights.
     *
     * @param n The number of assets.
     * @return A vector of equal asset weights.
     */
    public static double[] equalWeights (int n) {
        return DMatrixUtils.repeat(1.0 / n, n);
    }

    /**
     * Computes the volatility weighted asset weights for the given standard deviation vector.
     *
     * @param stddevs A vector of per-asset standard deviations.
     * @return A vector of volatility weighted asset weights.
     */
    public static double[] volatilityWeights(double[] stddevs) {
        // Initialize the return value:
        final double[] retval = new double[stddevs.length];

        // Compute the mean of standard deviations:
        final double meanStdDev = new DescriptiveStatistics(stddevs).getMean();

        // Populate the return value:
        for (int i = 0; i < retval.length; i++) {
            retval[i] = meanStdDev / stddevs[i];
        }

        // Compute the sum of values:
        final double sum = DMatrixUtils.sum(retval);

        // Create a real vector for convenience:
        final RealVector vec = new ArrayRealVector(retval, false);

        // Normalize and return:
        return vec.mapDivide(sum).toArray();
    }

    /**
     * Adds signs to each double in the vector by element-by-element multiplication.
     *
     * @param vector The vector to be signified.
     * @param signs Signs.
     * @return A signified vector.
     */
    public static double[] addSigns (double[] vector, double[] signs) {
        // TODO: Should we multiple by the signs of signs?
        return MatrixUtils.createRealVector(vector).ebeMultiply(new ArrayRealVector(signs, false)).toArray();
    }

    /**
     * Adds signs to elements of the matrix row-by-row.
     *
     * @param matrix The matrix to be signified.
     * @param signs Signs.
     * @return A signigied matrix.
     */
    public static RealMatrix addSigns(RealMatrix matrix, double[] signs) {
        // Initialize the return value:
        double[][] retval = new double[matrix.getRowDimension()][];

        // Iterate over rows and add signs:
        for (int i = 0; i < retval.length; i++) {
            retval[i] = FinanceUtils.addSigns(matrix.getRow(i), signs);
        }

        // Done, return as a matrix:
        return MatrixUtils.createRealMatrix(retval);
    }

    /**
     * Returns signs of the input matrix.
     *
     * @param matrix The input matrix.
     * @return A matrix of signs.
     */
    public static double[][] getDirections (double[][] matrix) {
        // Create the return value:
        double[][] retval = new double[matrix.length][matrix[0].length];

        // Iterate over elements and populate the return value by signs:
        for (int i = 0; i < retval.length; i++) {
            for (int j = 0; j < retval[i].length; j++) {
                retval[i][j] = Math.signum(matrix[i][j]);
            }
        }

        // Done, return:
        return retval;
    }

    /**
     * Simulates directions with the given hit ratio.
     *
     * @param directions Actual directions.
     * @param hitRatio The hit-ratio.
     * @return A simulated directions matrix.
     */
    public static double[][] simulateDirections(double[][] directions, double hitRatio) {
        // Initialize the return value:
        double[][] retval = new double[directions.length][directions[0].length];

        // Iterate over each element and simulate:
        for (int i = 0; i < retval.length; i++) {
            for (int j = 0; j < retval[i].length; j++) {
                retval[i][j] = (Math.random() < hitRatio) ? directions[i][j] : -directions[i][j];
            }
        }

        // Done, return:
        return retval;
    }

    /**
     * Computes the portfolio returns with the given trade directions.
     *
     * @param returns The return series matrix.
     * @param weights The trade direction.
     * @return
     */
    public static double[] portfolioReturns(double[][] returns, double[][] weights) {
        // Initialize the return value:
        double[] retval = new double[weights.length];

        // Iterate over rows and compute per period returns:
        for (int row = 0; row < retval.length; row++) {
            retval[row] = FinanceUtils.computeReturns(returns[row], weights[row]);
        }

        // Done, return:
        return retval;
    }

    /**
     * Creates a sequence of exponentials.
     *
     * @param length The desired length.
     * @param lambda The lambda factor.
     * @return The series of exponentials.
     */
    public static double[] expSeq (int length, double lambda) {
        // Initialize the return value:
        double[] retval = new double[length];

        // Iterate along the length, compute and assign:
        for (int i = length - 1; i >= 0; i--) {
            // TODO: Could the pow be any faster?
            retval[length - i - 1] = Math.pow(lambda, i) * (1 - lambda);
        }

        // Done, return:
        return retval;
    }

    /**
     * Backfills a series.
     *
     * @param series The series to be backfilled.
     * @return A backfilled series.
     */
    public static double[][] backfillSeries(double[][] series) {
        // Check series:
        if (series.length == 0) {
            return series;
        }

        // Get a clone as a matrix:
        // TODO: In fact, no need to use matrix class. It will make things slower. Using now just for convenience (cloning and accurate access).
        RealMatrix matrix = MatrixUtils.createRealMatrix(series);

        // Iterate over columns and backfill:
        for (int col = 0; col < matrix.getColumnDimension(); col++) {
            for (int row = 1; row < matrix.getRowDimension(); row++) {
                if (matrix.getEntry(row, col) == Double.NaN) {
                    matrix.setEntry(row, col, matrix.getEntry(row - 1, col));
                }
            }
        }

        // Done, return the date:
        return matrix.getData();
    }
}
