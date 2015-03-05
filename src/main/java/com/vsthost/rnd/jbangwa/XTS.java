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

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Created by vst on 3/3/15.
 */
public class XTS {
    private final Calendar[] index;
    private final String[] headers;
    private final RealMatrix data;
    private final Map<Calendar, Integer> lookupIndex;
    private final Map<String, Integer> lookupHeader;
    private final Calendar[] sortedIndex;

    public static final SimpleDateFormat DATETIME_FORMAT = new SimpleDateFormat("yyyy-MM-dd");

    public XTS(Calendar[] index, String[] headers, RealMatrix data) {
        // Check inputs:
        if (index.length != data.getRowDimension() || headers.length != data.getColumnDimension()) {
            throw new RuntimeException(
                String.format("Dimensions do not match. Indices: %d, Headers: %d, Matrix: %dx%d", index.length, headers.length, data.getRowDimension(), data.getColumnDimension())
            );
        }

        // Save stuff:
        this.index = index;
        this.headers = headers;
        this.data = data;

        // Create a map for fast index lookup:
        this.lookupIndex = new HashMap<>();
        for (int i = 0; i < this.index.length; i++) {
            this.lookupIndex.put(this.index[i], i);
        }

        // Create a map for fast header lookup:
        this.lookupHeader = new HashMap<>();
        for (int i = 0; i < this.headers.length; i++) {
            this.lookupHeader.put(this.headers[i], i);
        }

        // Create sorted index:
        this.sortedIndex = Arrays.copyOf(this.index, this.index.length);
        Arrays.sort(this.sortedIndex);
    }

    @Override
    public String toString() {
        String retval = "";
        for (int i = 0; i < this.index.length; i++) {
            retval += String.format("[%s] %s\n", DATETIME_FORMAT.format(this.index[i].getTime()), Arrays.toString(this.data.getRow(i)));
        }
        return retval;
    }

    public Calendar[] getIndex() {
        return index;
    }

    public String[] getIndexStrings() {
        // Declare and initialize the return value:
        String[] retval = new String[this.index.length];

        // Iterate and convert:
        for (int i = 0; i < this.index.length; i++) {
            retval[i] = DATETIME_FORMAT.format(this.index[i].getTime());
        }

        // Done, return:
        return retval;
    }

    public String[] getHeaders() {
        return headers;
    }

    public RealMatrix getData() {
        return data;
    }

    public int getLength () {
        return this.index.length;
    }

    public int getIndexOf (Calendar index) {
        return this.lookupIndex.get(index);
    }

    public int getIndexOf (String index) {
        return this.getIndexOf(XTS.parseCalendar(index));
    }

    public int getIndexOfFuzzy (Calendar index) {
        int i = 0;
        for (i = 0; i < this.sortedIndex.length; i++) {
            int comp = this.sortedIndex[i].compareTo(index);
            if (comp >= 0) {
                break;
            }
        }
        return i;
    }

    public int getIndexOfFuzzy (String index) {
        return this.getIndexOfFuzzy(XTS.parseCalendar(index));
    }

    public XTS getReturns () {
        // Compute returns:
        final RealMatrix returns = FinanceUtils.getReturns(this.data);

        // Done, return:
        return new XTS(Arrays.copyOfRange(this.index, 1, this.index.length), this.headers, returns);
    }

    public RealVector getRow (int index) {
        return this.data.getRowVector(index);
    }

    public RealVector getRow (Calendar index) {
        return this.getRow(this.lookupIndex.get(index));
    }

    public RealVector getRow (String index) {
        return this.getRow(XTS.parseCalendar(index));
    }

    public RealVector getColumn (int index) {
        return this.data.getRowVector(index);
    }

    public RealVector getColumn (String index) {
        return this.getColumn(this.lookupHeader.get(index));
    }

    public RealMatrix getColumns (int[] indices) {
        return this.data.getSubMatrix(IntStream.range(0, this.index.length).toArray(), indices);
    }

    public RealMatrix getColumns (String[] indices) {
        return this.getColumns(Arrays.stream(indices).mapToInt(e -> this.lookupHeader.get(e)).toArray());
    }

    public static XTS read (Reader reader) throws IOException {
        // Initialize the index:
        List<Calendar> index = new ArrayList<>();

        // Initialize headers:
        List<String> headers = new ArrayList<>();

        // Initialize the data:
        List<double[]> data = new ArrayList();

        // Parse and get the iterarable:
        Iterable<CSVRecord> records = CSVFormat.EXCEL.parse(reader);

        // Iterate over the records and populate return value:
        boolean flag = false;
        for (CSVRecord record : records) {
            if (!flag) {
                // Read in the headers:
                for (int i = 1; i < record.size(); i++) {
                    headers.add(record.get(i));
                }

                // Set the flag:
                flag = true;
            }
            else {
                // Create data row:
                double[] row = new double[record.size() - 1];

                // Iterate over records:
                for (int i = 0; i < record.size(); i++) {
                    if (i == 0) {
                        index.add(XTS.parseCalendar(record.get(i)));
                    }
                    else {
                        row[i - 1] = Double.parseDouble(record.get(i));
                    }
                }

                // Add to the raw data:
                data.add(row);
            }
        }

        // Convert the index to an array:
        Calendar[] indexArray = new Calendar[index.size()];
        index.toArray(indexArray);

        // Convert the index to an array:
        String[] headersArray = new String[headers.size()];
        headers.toArray(headersArray);

        // Convert the data to the matrix representation:
        double[][] dataMatrix = new double[data.size()][];
        data.toArray(dataMatrix);

        // Done, return the XTS instance:
        return new XTS(indexArray, headersArray, MatrixUtils.createRealMatrix(dataMatrix));
    }

    public static XTS read (String filepath) throws IOException {
        // Initilize the reader:
        Reader in = new FileReader(filepath);

        // Done, return:
        return XTS.read(in);
    }

    public static Calendar parseCalendar (String value) {
        // Split the date value:
        final String[] split = value.split("-");

        // Create and return calendar instance:
        return new GregorianCalendar(Integer.parseInt(split[0]), Integer.parseInt(split[1]) - 1, Integer.parseInt(split[2]));
    }
}
