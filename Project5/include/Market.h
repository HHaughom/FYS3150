// Header file for Market classes

#pragma once

#include <string>

class Market {
  protected:
    const int N, bins;
    const double m0, dm, maxMoney;
    double* agents;
    double* dist;
    double* meanDist;
  public:
    Market(int N, double m0, double dm, double maxMoney);
    double* getAgents() {return this->agents;}
    double* getDist() {return this->dist;}
    double* getMeanDist() {return this->meanDist;}
    int getBins() {return this->bins;}
    int getN() {return this->N;}

    virtual int transactions(double diffLimit, double saving);
    virtual void simulateMarket(int runs, double saving, double diffLimit);
    void writeToFile(std::string filename);
};


class NeighboursMarket : public Market {
  public:
    NeighboursMarket(int N, double m0, double dm, double maxMoney)
      : Market(N, m0, dm, maxMoney) {}
    int transactions(double diffLimit, double saving, double alfa);
    void simulateMarket(int runs, double saving, double diffLimit, double alfa);
};


class FormerTransactionsMarket : public Market {
  private:
    int** formerTransactions;
  public:
    FormerTransactionsMarket(int N, double m0, double dm, double maxMoney);
    int transactions(double diffLimit, double saving, double alfa, double gamma);
    void simulateMarket(int runs, double saving, double diffLimit, double alfa, double gamma);
};
