/*
 * Copyright (c) 2013 - Idiap Research Institute - http://www.idiap.ch/
 * 
 * See file COPYING for the licence associated to this software
 * 
 * This file is a part of the TopicModels4Video software version: academic-1.0
*/

#include <vector>
#include <fstream>


namespace MatlabUtils {
    using namespace std;
    void readFileLines(const char* filename, vector<string> &out) {
        string line;
        ifstream infile(filename, ios_base::in);
        while (getline(infile, line, '\n')) {
            out.push_back(line);
        }
        infile.close();
    }
    void splitLine(string in, string sep, vector<string> &out) {
        int sepSize = sep.size();
        string::size_type prevPos = 0;
        string::size_type pos = in.find_first_of(sep);
        while (prevPos != string::npos) {
            if (pos == prevPos + sepSize) {
                prevPos = pos == string::npos ? pos : pos + sepSize;
                pos = in.find_first_of(sep, prevPos);
                continue;
            }
            string val = in.substr(prevPos, pos == string::npos ? string::npos : pos - prevPos).c_str();
            out.push_back(val);
            prevPos = pos == string::npos ? pos : pos + sepSize;
            pos = in.find_first_of(sep, prevPos);
        }

    }

    void read(const char* filename, double *&mat, int &w, int &h, bool dropPossibleLastEmptys = true) {
        const char* sep = " \t";
        vector<string> file;
        readFileLines(filename, file);

        string s = file[0];
        string::size_type pos = s.find_first_of(sep);
        w = 0;
        while (pos != string::npos) {
            pos = s.find_first_of(sep, pos + 1);
            w++;
            if (dropPossibleLastEmptys && pos == s.size()-1) w--;
        }
        w++;

        h = file.size();

        if (dropPossibleLastEmptys && file[h-1].size() == 0) h--;
        mat = new double[w*h];

        int y = 0;
        vector<string>::iterator it;
        for (it = file.begin(); it < file.end() && y<h; it++, y++) {
            s = *it;

            int x = 0;
            string::size_type prevPos = 0;
            pos = s.find_first_of(sep);
            while (prevPos != string::npos && x<w) {
                if (pos == prevPos + 1) {
                    prevPos = pos;
                    pos = s.find_first_of(sep, pos + 1);
                    continue;
                }
                double val = atof(s.substr(prevPos, pos == string::npos ? string::npos : pos - prevPos).c_str());
                mat[y*w + x] = val;
                prevPos = pos;
                pos = s.find_first_of(sep, pos + 1);
                x++;
            }
        }
    }
    void write(const char* filename, double *mat, int w, int h) {
        const char* sep = " ";
        ofstream out(filename);
        for (int y=0; y<h; y++) {
            int x = 0;
            for (; x < w-1; x++)
                out << mat[x + y*w] << sep;
            out << mat[x + y*w] << endl;
        }
        out.close();
    }
    void free(double *mat) {
        delete[] mat;
    }
}


