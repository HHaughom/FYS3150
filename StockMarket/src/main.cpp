// Main made as executable for project 5

#include <iostream>
#include "../include/Market.h"
#include "../include/Timer.h"
#include "../include/tools.h"

using std::cout;
using std::cin;
using std::endl;

int main() {
  // Initializing variables
  int id = -1;
  int N, eq, runs;
  double m0, lambda, dm, maxMoney, diffLimit, alfa, gamma;
  std::string filename;

  while(id != 0) {
    cout << "\n0. Exit" << endl;
    cout << "1. Test Market class" << endl;
    cout << "2. Write distribution data from Market class to file" << endl;
    cout << "3. Test NeighboursMarket class" << endl;
    cout << "4. Write distribution data from NeighboursMarket class to file" << endl;
    cout << "5. Test FormerTransactionsMarket class" << endl;
    cout << "6. Write distribution from FormerTransactionsMarket class to file" << endl;
    cout << "Choose id: ";
    cin >> id;

    switch(id) {
      default: {
        cout << "Bad usage, id = " << id << endl;
        break;
      }
      case 0: {
        // Exit program
        break;
      }
      case 1: {
        /*
        Test of the simulateMarket function
        */
        // Getting parameters
        cout << "\nInput parameters:\nAgents = ";
        cin >> N;
        cout << "Saving quotient, lambda = ";
        cin >> lambda;
        cout << "Initial money, m0 = ";
        cin >> m0;
        cout << "Max money = ";
        cin >> maxMoney;
        cout << "Stepsize in distribution, dm = ";
        cin >> dm;
        cout << "Difference limit for reaching equilibrium, diffLimit = ";
        cin >> diffLimit;

        cout << "\nRunning ..." << endl;

        // Starting timer
        Timer timer;

        // Initializing object and simulating market
        Market test(N, m0, dm, maxMoney);
        eq = test.transactions(diffLimit, lambda);

        /*
        // Writing out data in agents to console
        double* agents = test.getAgents();
        for(int i = 0; i < test.getBins(); i++) {
          cout << agents[i] << endl;
        }
        delete [] agents;
        */
        break;
      }
      case 2: {
        /*
        Case for simulating market and writing data to file
        */
        // Looping as long as user want to
        std::string ans = "Y";
        while(ans == "Y" || ans == "y") {
          // Getting parameters
          cout << "\nInput parameters:\nAgents = ";
          cin >> N;
          cout << "Saving quotient, lambda = ";
          cin >> lambda;
          cout << "Initial money, m0 = ";
          cin >> m0;
          cout << "Max money = ";
          cin >> maxMoney;
          cout << "Stepsize in distribution, dm = ";
          cin >> dm;
          cout << "Difference limit for reaching equilibrium, diffLimit = ";
          cin >> diffLimit;
          cout << "Runs = ";
          cin >> runs;

          cout << "\nRunning ..." << endl;

          // Initializing object and simulating market
          Market test(N, m0, dm, maxMoney);
          test.simulateMarket(runs, lambda, diffLimit);

          // Taking filename as input and writhig to file
          cout << "\nFilename of data: ";
          cin >> filename;
          test.writeToFile(filename);

          // Checking if user want to run again
          cout << "Do you want to run again? [Y/N]: ";
          cin >> ans;
        }
        break;
      }
      case 3: {
        /*
        Test of the simulateMarket function
        */
        // Getting parameters
        cout << "\nInput parameters:\nAgents = ";
        cin >> N;
        cout << "Saving quotient, lambda = ";
        cin >> lambda;
        cout << "Initial money, m0 = ";
        cin >> m0;
        cout << "Max money = ";
        cin >> maxMoney;
        cout << "Stepsize in distribution, dm = ";
        cin >> dm;
        cout << "Difference limit for reaching equilibrium, diffLimit = ";
        cin >> diffLimit;
        cout << "Alfa = ";
        cin >> alfa;

        cout << "\nRunning ..." << endl;

        // Starting timer
        Timer timer;

        // Initializing object and simulating market
        NeighboursMarket testNeighbours(N, m0, dm, maxMoney);
        eq = testNeighbours.transactions(diffLimit, lambda, alfa);


        // Writing out data in agents to console
        double* agents = testNeighbours.getAgents();
        for(int i = 0; i < testNeighbours.getBins(); i++) {
          cout << agents[i] << endl;
        }
        delete [] agents;
        break;
      }
      case 4: {
        // Testing NeighboursMarket
        // Looping as long as user want to
        std::string ans = "Y";
        while(ans == "Y" || ans == "y") {
          // Getting parameters
          cout << "\nInput parameters:\nAgents = ";
          cin >> N;
          cout << "Saving quotient, lambda = ";
          cin >> lambda;
          cout << "Initial money, m0 = ";
          cin >> m0;
          cout << "Max money = ";
          cin >> maxMoney;
          cout << "Stepsize in distribution, dm = ";
          cin >> dm;
          cout << "Difference limit for reaching equilibrium, diffLimit = ";
          cin >> diffLimit;
          cout << "Runs = ";
          cin >> runs;
          cout << "alfa = ";
          cin >> alfa;

          cout << "\nRunning ..." << endl;

          // Initializing object and simulating market
          NeighboursMarket neighbours(N, m0, dm, maxMoney);
          neighbours.simulateMarket(runs, lambda, diffLimit, alfa);

          // Taking filename as input and writhig to file
          cout << "\nFilename of data: ";
          cin >> filename;
          neighbours.writeToFile(filename);

          // Checking if user want to run again
          cout << "Do you want to run again? [Y/N]: ";
          cin >> ans;
        }
        break;
      }
      case 5: {
        /*
        Test of the simulateMarket function
        */
        // Getting parameters
        cout << "\nInput parameters:\nAgents = ";
        cin >> N;
        cout << "Saving quotient, lambda = ";
        cin >> lambda;
        cout << "Initial money, m0 = ";
        cin >> m0;
        cout << "Max money = ";
        cin >> maxMoney;
        cout << "Stepsize in distribution, dm = ";
        cin >> dm;
        cout << "Difference limit for reaching equilibrium, diffLimit = ";
        cin >> diffLimit;
        cout << "Alfa = ";
        cin >> alfa;
        cout << "Gamma = ";
        cin >> gamma;

        cout << "\nRunning ..." << endl;

        // Starting timer
        Timer timer;

        // Initializing object and simulating market
        FormerTransactionsMarket testFormerTransactions(N, m0, dm, maxMoney);
        eq = testFormerTransactions.transactions(diffLimit, lambda, alfa, gamma);


        // Writing out data in agents to console
        double* agents = testFormerTransactions.getAgents();
        for(int i = 0; i < testFormerTransactions.getBins(); i++) {
          cout << agents[i] << endl;
        }
        delete [] agents;
        break;
      }
      case 6: {
        // Testing FormerTransactionsMarket
        // Looping as long as user want to
        std::string ans = "Y";
        while(ans == "Y" || ans == "y") {
          // Getting parameters
          cout << "\nInput parameters:\nAgents = ";
          cin >> N;
          cout << "Saving quotient, lambda = ";
          cin >> lambda;
          cout << "Initial money, m0 = ";
          cin >> m0;
          cout << "Max money = ";
          cin >> maxMoney;
          cout << "Stepsize in distribution, dm = ";
          cin >> dm;
          cout << "Difference limit for reaching equilibrium, diffLimit = ";
          cin >> diffLimit;
          cout << "Runs = ";
          cin >> runs;
          cout << "alfa = ";
          cin >> alfa;
          cout << "gamma = ";
          cin >> gamma;

          cout << "\nRunning ..." << endl;

          // Initializing object and simulating market
          FormerTransactionsMarket former(N, m0, dm, maxMoney);
          former.simulateMarket(runs, lambda, diffLimit, alfa, gamma);

          // Taking filename as input and writhig to file
          cout << "\nFilename of data: ";
          cin >> filename;
          former.writeToFile(filename);

          // Checking if user want to run again
          cout << "Do you want to run again? [Y/N]: ";
          cin >> ans;
        }
        break;
      }
    }
  }
  return 0;
}
