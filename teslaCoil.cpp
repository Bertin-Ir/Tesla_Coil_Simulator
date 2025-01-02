#include<iostream>
#include<cmath>
#include<fstream>
#define PI acos(-1.0) 
#define coupling_coefficient 0.99 

using namespace std;

// Class representing a power supply
class PowerSupply {
    private:
        double maxVoltage; // Maximum voltage output
        double f; // Frequency of the voltage

    public:
        // Constructor to initialize the power supply parameters
        PowerSupply(double maxvoltage, double frequency) {
            maxVoltage = maxvoltage;
            f = frequency;
        }

        // Function to calculate voltage at a given time
        double calculateVoltage(double time) {
            return (maxVoltage * sin(2 * PI * f * time)); // Sine wave voltage calculation
        }
};

// Class representing a transformer
class Transformer {
    private:
        double primaryInductance; // Inductance of the primary coil
        double ratio; // Turns ratio of the transformer

    public:
        // Constructor to initialize transformer parameters
        Transformer(double primaryInductance, double ratio) {
            this->primaryInductance = primaryInductance;
            this->ratio = ratio;
        }

        // Function to step up the voltage from the power supply
        double stepUpVoltage(PowerSupply inputVoltage, double time) {
            return (inputVoltage.calculateVoltage(time) * ratio); // Transformed voltage
        }
};

// Class representing a capacitor in the Tesla Coil circuit
class Capacitor {
    private:
        double capacitance; // Capacitance value
        double voltage; // Voltage across the capacitor

    public:
        // Constructor to initialize capacitor parameters
        Capacitor(double capacitance, double voltage) {
            this->capacitance = capacitance;
            this->voltage = voltage;
        }

        // Function to get the current voltage of the capacitor
        double getVoltage() {
            return voltage; // Return current voltage
        }

        // Function to calculate the current flowing through the capacitor
        double calculateCurrent(Transformer vsource, PowerSupply inputv, double tStep, double time) {
            double voltagechange = voltage - vsource.stepUpVoltage(inputv, time); // Voltage difference
            return (capacitance * voltagechange / tStep); // Current calculation
        }

        // Update the voltage of the capacitor based on current
        void updateVoltage(double tStep, Transformer vsource, PowerSupply inputv, double time) {
            voltage += (tStep * calculateCurrent(vsource, inputv, tStep, time) / capacitance); // Update voltage
        }

        // Function to charge the capacitor using the Runge-Kutta method
        void charge(Transformer vsource, double seriesResistance, PowerSupply inputv, double tstep, double t) {
            // Calculate k values for the Runge-Kutta method
            double k1 = (vsource.stepUpVoltage(inputv, t) - voltage) / (seriesResistance * capacitance);
            double k2 = (vsource.stepUpVoltage(inputv, t) - (voltage + 0.5 * k1 * tstep)) / (seriesResistance * capacitance);
            double k3 = (vsource.stepUpVoltage(inputv, t) - (voltage + 0.5 * k2 * tstep)) / (seriesResistance * capacitance);
            double k4 = (vsource.stepUpVoltage(inputv, t) - (voltage + k3 * tstep)) / (seriesResistance * capacitance);
            voltage += tstep * (k1 + 2 * k2 + 2 * k3 + k4) / 6; // Update voltage using the average of k values
        }
};

// Class representing a spark gap in the circuit
class SparkGap {
    private:
        double onResistance; // Resistance when conducting
        double offResistance; // Resistance when not conducting
        double breakdownVoltage; // Voltage needed to initiate conduction
        double holdingVoltage; // Minimum voltage to keep conducting
        double arcRecoverytime; // Time it takes to recover after conducting
        double tlast; // Last time the spark gap conducted
        bool conducting; // State of the spark gap (conducting or not)

    public:
        // Constructor to initialize spark gap parameters
        SparkGap(double onR, double offR, double breakV, double hv, double arcT) {
            onResistance = onR;
            offResistance = offR;
            breakdownVoltage = breakV;
            holdingVoltage = hv;
            arcRecoverytime = arcT;
            tlast = 0.0; // Initialize last conduction time
            conducting = false; // Start as not conducting
        }

        // Function to get the resistance of the spark gap based on voltage
        double getResistance(Capacitor vc, double t) {
            // Check if the spark gap should conduct
            if (vc.getVoltage() >= breakdownVoltage && (t - tlast) > arcRecoverytime) {
                conducting = true; // Spark gap conducts
                tlast = t; // Update last conduction time
                return onResistance; // Return on resistance
            }
            // If not conducting, return off resistance
            else if (vc.getVoltage() < holdingVoltage && !conducting) {
                conducting = false; // Update state to not conducting
                return offResistance; // Return off resistance
            }
            return offResistance; // Default return to avoid undefined behavior
        }

        // Function to calculate the current across the spark gap
        double SparkCurrent(Capacitor vc, Transformer vtrans, PowerSupply inputv, double tstep, double t) {
            return (vc.calculateCurrent(vtrans, inputv, tstep, t) * (vc.getVoltage() - vtrans.stepUpVoltage(inputv, t)) / getResistance(vc, t)); // Current calculation
        }

        // Function to calculate voltage across the spark gap
        double SparkVoltage(Capacitor vc, Transformer vtrans, PowerSupply inputv, double tstep, double t) {
            return (SparkCurrent(vc, vtrans, inputv, tstep, t) * getResistance(vc, t)); // Voltage calculation
        }
};

// Class representing both the primary and secondary coils of the Tesla Coil
class Coil {
    private:
        double inductance; // Inductance of the coil
        double resistance; // Resistance of the coil
        bool isprimary; // Flag to indicate if it's a primary coil

    public:
        // Constructor to initialize coil parameters
        Coil(double inductance, double resistance, bool isprimary) {
            this->inductance = inductance;
            this->resistance = resistance;
            this->isprimary = isprimary; // Set primary flag
        }

        // Getters for inductance and resistance
        double getInductance() {
            return inductance; // Return inductance
        }
        double getResistance() {
            return resistance; // Return resistance
        }

        // Function to calculate voltage induced in the coil
        double calculateVoltage(SparkGap I_change, Capacitor vc, Transformer vtrans, 
            PowerSupply inputv, double tstep, double t, Coil secondary) {
            if (isprimary) {
                return (inductance * I_change.SparkCurrent(vc, vtrans, inputv, tstep, t) / tstep); // Primary voltage calculation
            } else {
                double M = coupling_coefficient * sqrt(inductance * secondary.getInductance()); // Mutual inductance calculation
                return (-M * I_change.SparkCurrent(vc, vtrans, inputv, tstep, t) / tstep); // Secondary voltage calculation
            }
        }
};

// Main class to simulate the entire Tesla Coil system
class TeslaCoil {
    private:
        PowerSupply* p; // Pointer to PowerSupply object
        Transformer* Trans; // Pointer to Transformer object
        SparkGap* spark; // Pointer to SparkGap object
        Capacitor* cap; // Pointer to Capacitor object
        Coil* PrimaryCoil, *SecondaryCoil; // Pointers to Coil objects

    public:
        // Constructor to initialize the Tesla Coil components
        TeslaCoil(PowerSupply *pow, Transformer* tran, SparkGap* spa, 
            Coil * primary, Coil* secondary, Capacitor* c) {
            p = pow; // Initialize power supply
            Trans = tran; // Initialize transformer
            spark = spa; // Initialize spark gap
            cap = c; // Initialize capacitor
            PrimaryCoil = primary; // Initialize primary coil
            SecondaryCoil = secondary; // Initialize secondary coil
        }

        // Method to simulate the Tesla Coil operation
        void simulate(double dt, double simulationTime, ofstream& allData, ofstream& voltTimeData) {
            double time = 0; // Current simulation time

            // Write headers for the data files
            allData << "time," << "Input Voltage," << "Transformed Voltage," << "Capacitor Voltage,"
                    << "Spark Resistance," << "Primary coil Voltage," << "Secondary Coil voltage,"
                    << "Capacitor Current," << "Spark Current," << "Spark Voltage" << endl;

            voltTimeData << "time," << "Input Voltage," << "Transformed Voltage," << "Capacitor Voltage,"
                         << "Primary coil Voltage," << "Secondary Coil Voltage," << "Spark Voltage" << endl;

            // Simulation loop
            while (time <= simulationTime) {
                double inputVoltage = p->calculateVoltage(time); // Calculate input voltage
                double transformedVoltage = Trans->stepUpVoltage(*p, time); // Calculate transformed voltage
                double capCurrent = cap->calculateCurrent(*Trans, *p, dt, time); // Calculate capacitor current
                cap->updateVoltage(dt, *Trans, *p, time); // Update capacitor voltage
                cap->charge(*Trans, 100, *p, dt, time); // Charge capacitor
                double capacitorVoltage = cap->getVoltage(); // Get updated capacitor voltage
                double sparkResistance = spark->getResistance(*cap, time); // Get spark gap resistance
                double sparkCurr = spark->SparkCurrent(*cap, *Trans, *p, dt, time); // Calculate spark current
                double sparkV = spark->SparkVoltage(*cap, *Trans, *p, dt, time); // Calculate spark voltage
                double primaryVoltage = PrimaryCoil->calculateVoltage(*spark, *cap, *Trans, 
                                                                      *p, dt, time, *SecondaryCoil); // Calculate primary coil voltage
                double secondaryVoltage = SecondaryCoil->calculateVoltage(*spark, *cap, *Trans, 
                                                                          *p, dt, time, *PrimaryCoil); // Calculate secondary coil voltage

                // Save all parameters to allData file
                allData << time << "," << inputVoltage << "," << transformedVoltage << "," 
                        << capacitorVoltage << "," << sparkResistance << "," 
                        << primaryVoltage << "," << secondaryVoltage << "," 
                        << capCurrent << "," << sparkCurr << "," << sparkV << endl;

                // Save voltage and time data to voltTimeData file
                voltTimeData << time << "," << inputVoltage << "," << transformedVoltage << "," << capacitorVoltage << ","
                             << primaryVoltage << "," << secondaryVoltage << "," << sparkV << endl;

                time += dt; // Increment time
            }
        }

        // Function to display simulation results
        void displaySimulationResults() {
            cout << "Simulation Results: \n";
            ifstream alldata; // Input file stream for results
            alldata.open("all_components_data.csv"); // Open data file
            if (alldata.fail()) {
                cerr << "Failed to open the file\n"; // Error message if file fails to open
                exit(1);
            }
            char ch;

            // Read and display the contents of the data file
            while (!alldata.eof()) {
                alldata.get(ch);
                if (!alldata.eof())
                    cout << ch; // Output each character
            }

            alldata.close(); // Close the file
        }
};

int main() 
{ 
    // Define the simulation parameters 
    double dt; // Time step for the simulation 
    double simulationTime; // Total time to simulate 

    cout << "Enter the simulation duration: "; 
    cin >> simulationTime; // Read simulation duration

    cout << "Enter the simulation step size: "; 
    cin >> dt; // Read time step size

    // Define the component parameters 
    double powerSupplyMaxVoltage = 120.0; // Max Voltage for power supply 
    double powerSupplyFrequency = 60.0;   // Frequency for power supply 

    double transformerPrimaryInductance = 10.0; // Inductance of the transformer primary coil in Henries 
    double transformerRatio = 500.0;            // Turns ratio of the transformer 
 
    double sparkGapOnResistance = 1.0;           // Resistance when spark gap is on (Ohms) 
    double sparkGapOffResistance = 1e9;          // Resistance when spark gap is off (Gigaohms) 
    double sparkGapBreakdownVoltage = 1000.0;    // Voltage to breakdown the spark gap (Volts) 
    double sparkGapHoldingVoltage = 0.001;       // Minimum voltage to keep spark gap conducting (Volts) 
    double sparkGapArcRecoveryTime = 0.0001;     // Arc recovery time (Seconds) 

    double capacitorCapacitance = 20e-6; // Capacitance value (20 uF) 
    double capacitorInitialVoltage = 1.0; // Initial voltage across capacitor (Volts) 

    double primaryCoilInductance = 31.6628e-6; // Inductance of primary coil (31.6628 uH) 
    double primaryCoilResistance = 100e-3;     // Resistance of primary coil (100 mOhm) 

    double secondaryCoilInductance = 1e-3; // Inductance of secondary coil (1 mH) 
    double secondaryCoilResistance = 100e-3; // Resistance of secondary coil (100 mOhm) 

    // Create the components for the Tesla Coil
    PowerSupply powerSupply(powerSupplyMaxVoltage, powerSupplyFrequency); 
    Transformer stepUpTransformer(transformerPrimaryInductance, transformerRatio); 
    SparkGap sparkGap(sparkGapOnResistance, sparkGapOffResistance, 
                      sparkGapBreakdownVoltage, sparkGapHoldingVoltage, sparkGapArcRecoveryTime); 
    Capacitor primaryCapacitor(capacitorCapacitance, capacitorInitialVoltage); 
    Coil primaryCoil(primaryCoilInductance, primaryCoilResistance, true); 
    Coil secondaryCoil(secondaryCoilInductance, secondaryCoilResistance, false); 

    // Create the TeslaCoil object, passing references to the components 
    TeslaCoil teslaCoil(&powerSupply, &stepUpTransformer, &sparkGap, 
                         &primaryCoil, &secondaryCoil, &primaryCapacitor); 

    // Declare output files 
    ofstream allData("all_components_data.csv"); 
    ofstream voltTimeData("voltage_time.csv"); 

    // Check if output files are open 
    if (!allData.is_open()) { 
        std::cerr << "Failed to open all_components_data.csv" << std::endl; 
        return EXIT_FAILURE; 
    } 

    if (!voltTimeData.is_open()) { 
        std::cerr << "Failed to open voltage_time.csv" << std::endl; 
        return EXIT_FAILURE; 
    } 

    // Run the simulation 
    teslaCoil.simulate(dt, simulationTime, allData, voltTimeData); 

    // Close the output files 
    allData.close(); 
    voltTimeData.close(); 

    // Display the simulation results 
    teslaCoil.displaySimulationResults(); 

    return 0; 
}