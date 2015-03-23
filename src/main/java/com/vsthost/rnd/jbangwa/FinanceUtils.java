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

import com.vsthost.rnd.commons.math.ext.linear.EMatrixUtils;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.Arrays;

/**
 * Created by vst on 27/2/15.
 */
public class FinanceUtils {

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

    public static RealMatrix getReturns (final RealMatrix prices) {
        return org.apache.commons.math3.linear.MatrixUtils.createRealMatrix(FinanceUtils.getReturns(prices.getData()));
    }


    public static RealMatrix detrendReturns(final RealMatrix returns) {
        // Subtract col means from rows and return:
        return EMatrixUtils.rowSubtract(returns, EMatrixUtils.colMeans(returns));
    }

    public static double computeReturns (double[] returns, double[] weights) {
        double total = 0.0;
        for (int i = 0; i < returns.length; i++) {
            total += returns[i] * weights[i];
        }
        return total;
    }

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

    public static double[] portfolioReturns(final RealMatrix returns) {
        // Compute equal weights and call overload:
        return FinanceUtils.portfolioReturns(returns, FinanceUtils.equalWeights(returns.getColumnDimension()));
    }

    public static double[] equalWeights (int n) {
        double[] weights = new double[n];
        Arrays.fill(weights, 1.0 / n);
        return weights;
    }

    public static double[] volatilityWeights(double[] doubles) {
        double[] retval = new double[doubles.length];
        double meanStdDev = new DescriptiveStatistics(doubles).getMean();
        for (int i = 0; i < retval.length; i++) {
            retval[i] = meanStdDev / doubles[i];
        }
        RealVector vec = MatrixUtils.createRealVector(retval);
        return vec.mapDivide(new DescriptiveStatistics(retval).getSum()).toArray();
    }

    public static double[] addSigns (double[] vector, double[] signs) {
        return MatrixUtils.createRealVector(vector).ebeMultiply(MatrixUtils.createRealVector(signs)).toArray();
    }

    public static RealMatrix addSigns(RealMatrix matrix, double[] signs) {
        double[][] retval = new double[matrix.getRowDimension()][];
        for (int i = 0; i < retval.length; i++) {
            retval[i] = FinanceUtils.addSigns(matrix.getRow(i), signs);
        }
        return MatrixUtils.createRealMatrix(retval);
    }

    public static double[][] getDirections (double[][] returns) {
        double[][] retval = new double[returns.length][returns[0].length];
        for (int i = 0; i < retval.length; i++) {
            for (int j = 0; j < retval[i].length; j++) {
                retval[i][j] = Math.signum(returns[i][j]);
            }
        }
        return retval;
    }

    public static double[][] simulateDirections(double[][] directions, double hitRatio) {
        double[][] retval = new double[directions.length][directions[0].length];
        for (int i = 0; i < retval.length; i++) {
            for (int j = 0; j < retval[i].length; j++) {
                if (Math.random() < hitRatio) {
                    retval[i][j] = directions[i][j];
                }
                else {
                    retval[i][j] = -directions[i][j];
                }
            }
        }
        return retval;
    }

    public static double[] portfolioReturns(double[][] returnsOS, double[][] simulatedPortfolioDirections) {
        double[] retval = new double[simulatedPortfolioDirections.length];
        for (int row = 0; row < retval.length; row++) {
            retval[row] = FinanceUtils.computeReturns(returnsOS[row], simulatedPortfolioDirections[row]);
        }
        return retval;
    }

    public static double[] expSeq (int length, double lambda) {
        double[] retval = new double[length];
        for (int i = length - 1; i >= 0; i--) {
            retval[length - i - 1] = Math.pow(lambda, i) * (1 - lambda);
        }
        return retval;
    }

    public static double[][] backfillSeries(double[][] prices) {
        // Check prices:
        if (prices.length == 0) {
            return prices;
        }

        // Get a clone as a matrix:
        RealMatrix matrix = MatrixUtils.createRealMatrix(prices);

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
