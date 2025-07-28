
#define FBA_INFO_FILE "${FBA_INFO_FILE}"

/**FBAProcessor Class Definition**/

constexpr double INTEGRATION_STEP_H = 1;   // Δt 


// Utility: controlla se la stringa s termina con il suffisso suffix
static bool endsWith(const std::string& s, const std::string& suffix) {
    if (suffix.size() > s.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

/**
 * @class FBAProcessor
 * @brief A singleton class that manages flux balance analysis processes.
 *
 * The FBAProcessor class is designed using the Singleton design pattern to ensure
 * that only one instance of this class is created throughout the lifetime of the application.
 * This class handles the initialization, management, and computation of flux balance analysis
 * based on metabolic reaction data.
 */
class FBAProcessor {
		public:
    /**
     * @brief Retrieves the single instance of the FBAProcessor class.
     * 
     * This method ensures that the FBAProcessor class is instantiated only once
     * and the same instance is returned with every call. The instance is guaranteed
     * to be destroyed only when the application exits, ensuring controlled and predictable
     * lifecycle management of the singleton resource.
     * 
     * @return Reference to the singleton FBAProcessor instance.
     */
    static FBAProcessor& getInstance() {
        static FBAProcessor instance; // Static variable ensures that the instance is created only once.
        return instance;
    }

    /**
     * @brief Deletes the copy constructor to prevent copying of the singleton instance.
     */
    FBAProcessor(FBAProcessor const&) = delete;
    /**
     * @brief Deletes the assignment operator to prevent assignment and copying of the singleton instance.
     */
    void operator=(FBAProcessor const&) = delete;



    /**
     * @brief Processes metabolic changes and calculates the rate of a specified transition.
     * 
     * This method serves as the primary access point for external components to interact
     * with the FBAProcessor. It ensures that the system is initialized before proceeding
     * with updates and calculations. The method handles the following tasks:
     * - Initializes the data structures if not already initialized.
     * - Updates the bounds based on current metabolic values and solves the LP problems.
     * - Computes the rate of the specified transition based on the new LP solution.
     * 
     * @param Value Pointer to the array containing current metabolic values.
     * @param vec_fluxb Reference to a vector of LPprob objects representing flux balance problems.
     * @param NumTrans Mapping of transition names to their corresponding indices.
     * @param NumPlaces Mapping of place names to their corresponding indices in the metabolic array.
     * @param NameTrans Vector containing names of transitions, indexed by transition IDs.
     * @param Trans Pointer to the structure containing transition information (not used directly in this method but may be used for extending functionality).
     * @param T The index of the transition for which the rate is to be calculated.
     * @param time The current simulation time, used to check if updates are needed based on time-driven changes.
     * @return The rate of the specified transition after processing the metabolic changes.
     * This rate could be used to determine the speed or likelihood of the transition occurring.
     */
		double process(double *Value,
                             vector<FBGLPK::LPprob>& vec_fluxb,
                             map<string,int>& NumTrans,
                             map<string,int>& NumPlaces,
                             const vector<string>& NameTrans,
                             const InfTr* Trans,
                             int T,
                             const double& time)
		{
				// --- init on first call ---
				if (!init) {
				    init_data_structures_class(NameTrans, vec_fluxb, Value, NumPlaces);
				    mapReactionsFromProblems(vec_fluxb);
				    mapMetabolitesToReactions();
				    loadAndApplyFBAReactionUpperBounds(vec_fluxb, "EX_upper_bounds_FBA.csv");
				    updateNonFBAReactionUpperBoundsFromFile(vec_fluxb, "EX_upper_bounds_nonFBA.csv");
				    loadMuMaxValues(vec_fluxb);
				    loadGeneRules(vec_fluxb, "GeneRules.txt");
				    debugPrintGeneRules();
				    init = true;
				    std::cout << "[DEBUG process] Initialized FBAProcessor\n";
				}

				// --- aggiorno deadBacterialSpecies ---
				for (const auto& kv : FBAproblems) {
				    size_t idx = kv.second;
				    if (hasMultiSpecies) {
				        double pop = std::floor(Value[ NumPlaces.at(problemBacteriaPlace.at(idx)) ]);
				        if (pop < 1) {
				            deadBacterialSpecies[idx] = true;
				        }
				    }
				}

					bool newStep = std::isnan(stepTime) || (std::fabs(time - stepTime) > tol);
					if (newStep && !std::isnan(lastSolvedTime) && std::fabs(time - lastSolvedTime) < minDtFBA)
							newStep = false;                              

					if (newStep) {

							std::cout << "[DEBUG newStep] ▶ Entering FBA update block at time = " << time
										    << " (previous stepTime = " << stepTime << ")\n";

							stepTime       = time;
							lastSolvedTime = time; 

							auto changed = checkSignificantChange(Value, NumPlaces, time);
							std::cout << "[DEBUG newStep]   changed metabolites (" << changed.size() << "): ";
							for (const auto& m : changed) std::cout << m << ", ";
							std::cout << "\n";

							reactionsToUpdate.clear();
							for (auto& met : changed) {
									auto it = metaboliteToReactions.find(met);
									if (it != metaboliteToReactions.end()) {
										  for (auto& rxn : it->second)
										      reactionsToUpdate.insert(rxn);
									}
							}
							std::cout << "[DEBUG newStep]   reactionsToUpdate size = "
										    << reactionsToUpdate.size() << "\n";

							bool bioUpdated = false;
							for (auto& bioPlace : problemBiomassPlace) {
									if (changed.count(bioPlace.second)) {
										  std::cout << "[DEBUG newStep]   Biomass place '"
										            << bioPlace.second << "' changed → updating biomass bounds\n";
										  updateAllBiomassReactionsUpperBounds(Value, NumPlaces, vec_fluxb);
										  bioUpdated = true;
										  break;
									}
							}
							if (!bioUpdated) {
									std::cout << "[DEBUG newStep]   No biomass update needed\n";
							}

							std::cout << "[DEBUG newStep]   Calling updateFluxBoundsAndSolve …\n";
							updateFluxBoundsAndSolve(Value, vec_fluxb, NumPlaces, time);
							std::cout << "[DEBUG newStep]   Returned from updateFluxBoundsAndSolve\n";

							updateConcentrations(Value, NumPlaces, changed);
							std::cout << "[DEBUG newStep]   Previous concentrations updated\n";

							std::cout << "[DEBUG newStep] ◀ Exiting FBA update block\n\n";
					} else {
							std::cout << "[DEBUG process] Skipping FBA updates "
										    << "(newStep=" << std::boolalpha << newStep
										    << ", stepTime=" << stepTime
										    << ", lastSolvedTime=" << lastSolvedTime << ")\n";
					}


				// --- calcolo il rate per la sola transizione T ---
				size_t problemIndex = FBAproblems[ NameTrans[T] ];
				if (deadBacterialSpecies.count(problemIndex)) {
				    std::cout << "[DEBUG process] Species dead → rate=0\n";
				    return 0.0;
				}
				if (transitionsWithoutFile.count(NameTrans[T])) {
				    std::cout << "[DEBUG process] No-file transition → rate=0\n";
				    return 0.0;
				}

				std::cout << "[DEBUG process] computeRate for '" 
				          << NameTrans[T] << "'\n";
				double rate = computeRate(vec_fluxb, NumPlaces, NameTrans, Value, T, decimalTrunc, time);
				std::cout << "[DEBUG process] rate = " << rate << "\n";
				return rate;
		}




		private:
		
		/* -------------- Simulation (ODE + FBA) parameters (used for efficent update of fba bounds) -------------- */
		double stepTime       = std::numeric_limits<double>::quiet_NaN(); 
		double lastSolvedTime = std::numeric_limits<double>::quiet_NaN(); 
		const  double tol     = 1e-12;     
		const  double minDtFBA = 0.0;       
		
		double lastFbaTime = std::numeric_limits<double>::lowest();
		// place → set di reactionName che lo coinvolgono
		std::unordered_map<std::string, std::set<std::string>> metaboliteToReactions;
		// species‑specific reaction → base reaction
		unordered_map<string,string> speciesToBaseRxn;
		
		unordered_set<std::string> projectedFBA;

    struct GeneRule {
        bool timeSpecified;       // true if a time condition is specified
        double time;              // simulation time condition (only valid if timeSpecified is true)

        bool placeSpecified;      // true if a metabolite condition is specified
        std::string place;        // name of the metabolite (e.g., "glc_D_e")

        bool thresholdSpecified;  // true if a threshold condition is specified
        double threshold;         // concentration threshold (only valid if thresholdSpecified is true)

        std::string compType;     // comparison type: "min", "max", or "NA"

        std::string reactionID;   // reaction identifier (e.g., "LACZ")
        double newLB;             // new lower bound for the reaction
        double newUB;             // new upper bound for the reaction
        std::string bacterium;    // bacterium name for which the rule is associated
        size_t lpIndex;           // index in the vector of LP problems (default: numeric_limits<size_t>::max())
        bool applied;             // flag to track if the rule has been applied
    };

    struct PairHash {
        size_t operator()(const std::pair<std::string, size_t>& p) const {
            // Combino gli hash
            // (puoi usare qualunque formula di combinazione)
            auto h1 = std::hash<std::string>()(p.first);
            auto h2 = std::hash<size_t>()(p.second);
            // Semplice combiner stile boost::hash_combine
            h1 ^= (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
            return h1;
        }
    };

    // Se vuoi un comparator personalizzato (ma di solito la pair ha già operator==)
    struct PairEq {
        bool operator()(const std::pair<std::string, size_t>& a,
                        const std::pair<std::string, size_t>& b) const {
            return (a.first == b.first && a.second == b.second);
        }
    };
    
		 // generic hash for <string,size_t>
		struct RxnLPHash {
				std::size_t operator()(const std::pair<std::string,std::size_t>& p) const noexcept {
				    return std::hash<std::string>()(p.first) ^ (p.second * 0x9e3779b97f4a7c15ULL);
				}
		};

     set<string> reactionsToUpdate;
     string firstTransitionName = ""; // Tracks problems associated with dead bacterial species
     unordered_map<size_t, bool> deadBacterialSpecies; // Tracks problems associated with dead bacterial species
    /**
     * comment.
     */
     unordered_map<size_t, bool> problemsWithLowBiomass; // Tracks problems with biomass below the minimum threshold
			const double Lcutoff = 1e-6; // Define a small cutoff value for limit the biomass upper bound
    /**
     * Unordered set to store transitions identified as related to biomass.
     */
     unordered_set<string> biomassTransitions;
			unordered_set<string> transitionsWithoutFile;
    std::unordered_map<
        std::pair<std::string, size_t>,
        double,
        PairHash,
        PairEq
    > NonFBAReactionBaseUB;
    
		std::unordered_map<
				std::pair<std::string,std::size_t>,
				std::string,
				PairHash,
				PairEq
		> ReactionSubtype;         // "exchange" | "internal"
		
	
    unordered_map<string, unordered_set<size_t>> reactionToFileMap; // Maps each reaction to a set of files where it appears

    unordered_map<int, double> bacteriaToBioMin; // Maps each reaction to a set of files where it appears

    unordered_map<int, double> bacteriaToBioMax; // Maps each reaction to a set of files where it appears

    unordered_map<int,double> bacteriaToBioMean; // Maps each reaction to a set of files where it appears

    /**
     * Multiplicative constant for the reaction.
     */
    double multiplicativeConstant = 1;

    /**
     * Comment.
     */
    bool hasMultiSpecies = false;

    /**
     * Comment.
     */
    double hasBioMASS = false;

    /**
     * Maps each transition to its corresponding reaction. Used for linking network transitions to specific biochemical reactions.
     */
    unordered_map<string, string> FBAreact;
    /**
     * Stores mappings from reactions to their places. This helps in maintaining the state of metabolic concentrations involved in reactions.
     */
    unordered_map<string, string> FBAplace;
    /**
     * Maps each reaction to an index of the corresponding LP problem in a vector. This facilitates quick access to the problem related to a particular reaction.
     */
    unordered_map<string, size_t> FBAproblems;
    /**
     * Maps each LP problem index to a set of reactions associated with that problem.
     * This helps in managing and optimizing the parallel processing of flux balance analysis problems,
     * ensuring that all reactions associated with a particular problem are processed together.
     */
    unordered_map<size_t, set<string>> problemsToReactions;
    /**
     * A set of all reactions that are part of the FBA model. Helps in ensuring that each reaction is processed uniquely.
     */
    set<string> FBAreactions;
    /**
     * A set of all problem index that are part of the FBA model.
     */
    set<size_t> problems;
    /**
     * Initialization status flag. Used to ensure that the FBAProcessor is set up before any operations are performed.
     */
    bool init;
    /**
     * Precision setting for numerical operations, specifically truncating numbers to a fixed number of decimal places.
     */
    const double decimalTrunc = 16;
    /**
     * Stores the last time the FBA model was updated, used to prevent unnecessary recalculations within the same time frame.
     */
     double FBAtime = -1;
    /**
     * Molecular Weight scaling factor.
     */
    double Mw = 1;
    double minDeltaThreshold = 1e-16;
    unordered_map<string, double> ReactionMultiplicity;
    /**
     * Pointer to an array of pointers, each pointing to a set of variables representing the results of the linear programming problems.
     */
    double** Vars;
    /**
     * Tracks the previous concentrations for all places, allowing for comparison against new values to detect significant changes.
     */
    unordered_map<string, double> previousConcentrations;
    /**
     * Comment.
     */
    unordered_map<size_t, string> problemBacteriaPlace;
    /**
     * Comment.
     */
    unordered_map<size_t, string> problemBiomassPlace;
    /**
     * Comment.
     */
    unordered_map<string, string> reactionToBacteria;
    /**
     * Comment.
     */
    unordered_map<string, string> reactionToBacteriaBIOMASS;
    /**
     * Threshold for determining significant changes in metabolic concentrations, expressed as a absolute value.
     */
    double absEpsilon = g_absEpsilon; // Initially set to 0%, can be configured.
    /**
     * Threshold for determining significant changes in metabolic concentrations, expressed as a percentage.
     */
    double relEpsilon = g_relEpsilon; // Initially set to 0%, can be configured.
    double scalingFactor = 1e-12; // Initially set to 0%, can be configured.
    
		double systemExternalVolume = 1.0;   
		bool   externalVolumeLoaded   = false; 
    /**
     * Counter to track how many times the network state has undergone significant changes.
     */
    double count = 0;
    /**
     * Enumeration to distinguish between input and output transitions within the metabolic network.
     */
    enum TransitionType { Input, Output };
    /**
     * Private constructor to prevent instantiation outside of the getInstance() method, ensuring singleton behavior.
     */
    FBAProcessor() : init(false) {
    }

    /**
     * Maps each metabolite to a set of problem indices that are affected by changes in this metabolite.
     * This helps optimize the processing by updating only relevant LP problems when specific metabolite concentrations change.
     */
    unordered_map<string, set<size_t>> metaboliteToProblems;

		/*Max possible flux value**/
		std::unordered_map<std::size_t,double> muMaxMap;


    std::vector<GeneRule> geneRules;
    bool geneRulesLoaded;

    /**
     * @brief Loads gene regulation rules from a file with a default name.
     *        If the file is not present, the simulation continues without gene rules.
     * @param vec_fluxb The vector of LP problems (used to map bacterium names to LP indices).
     * @param filename The file name to load; default is "GeneRules.txt".
     */
    void loadGeneRules(const vector<class FBGLPK::LPprob>& vec_fluxb, const std::string &filename = "GeneRules.txt") {
        std::ifstream infile(filename.c_str());
        if (!infile) {
            std::cout << "Gene rules file '" << filename 
                      << "' not found. Continuing simulation without gene rules." << std::endl;
            geneRulesLoaded = false;
            return;
        }
        geneRules.clear();
        std::string line;
        while (std::getline(infile, line)) {
            if (line.empty() || line[0] == '#') continue;
            size_t commentPos = line.find('#');
            if (commentPos != std::string::npos) line = line.substr(0, commentPos);
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            std::string token;
            while (std::getline(iss, token, ',')) {
                token.erase(0, token.find_first_not_of(" \t\r\n"));
                token.erase(token.find_last_not_of(" \t\r\n") + 1);
                if (!token.empty()) tokens.push_back(token);
            }
            if (tokens.size() < 7) {
                std::cerr << "Skipping invalid gene rule: " << line << std::endl;
                continue;
            }
            GeneRule rule;
            if (tokens[0] == "NA") { rule.timeSpecified = false; } else { rule.timeSpecified = true; rule.time = std::stod(tokens[0]); }
            if (tokens[1] == "NA") { rule.placeSpecified = false; } else { rule.placeSpecified = true; rule.place = tokens[1]; }
            if (tokens[2] == "NA") { rule.thresholdSpecified = false; } else { rule.thresholdSpecified = true; rule.threshold = std::stod(tokens[2]); }
            // tokens[3] is reactionID
            rule.reactionID = tokens[3];
            // tokens[4] is newLB, tokens[5] is newUB
            rule.newLB = std::stod(tokens[4]);
            rule.newUB = std::stod(tokens[5]);
            // tokens[6] is bacterium name
            rule.bacterium = tokens[6];
            // If a comparison type is provided as an optional 7th token, use it; otherwise, set to "NA"
            if (tokens.size() >= 8) {
                rule.compType = tokens[7];
            } else {
                rule.compType = "NA";
            }
            // Determine the LP index for the bacterium; assume findLPIndex is defined
            rule.lpIndex = findLPIndex(vec_fluxb, rule.bacterium);
            if (rule.lpIndex == std::numeric_limits<size_t>::max()) {
                std::cerr << "Error: Bacterium '" << rule.bacterium << "' not found among LP problems." << std::endl;
            }
            // Set the applied flag to false initially
            rule.applied = false;
            geneRules.push_back(rule);
        }
        infile.close();
        geneRulesLoaded = true;
        std::cout << "Loaded " << geneRules.size() << " gene rule(s) from " << filename << std::endl;
    }
    void debugPrintGeneRules() {
        std::cout << "DEBUG: Loaded " << geneRules.size() << " gene rule(s):" << std::endl;
        for (size_t i = 0; i < geneRules.size(); ++i) {
            const GeneRule &rule = geneRules[i];
            std::cout << "GeneRule[" << i << "]: ";
            std::cout << "Time = ";
            if (rule.timeSpecified)
                std::cout << rule.time;
            else
                std::cout << "NA";
            std::cout << ", Place = ";
            if (rule.placeSpecified)
                std::cout << rule.place;
            else
                std::cout << "NA";
            std::cout << ", Threshold = ";
            if (rule.thresholdSpecified)
                std::cout << rule.threshold;
            else
                std::cout << "NA";
            std::cout << ", ReactionID = " << rule.reactionID;
            std::cout << ", newLB = " << rule.newLB;
            std::cout << ", newUB = " << rule.newUB;
            std::cout << ", Bacterium = " << rule.bacterium;
            std::cout << ", lpIndex = " << rule.lpIndex;
            std::cout << ", compType = " << rule.compType;
            std::cout << ", applied = " << (rule.applied ? "YES" : "NO") << std::endl;
        }
    }
    void applyGeneRegulationRules(double *Value, vector<class FBGLPK::LPprob>& vec_fluxb, map<string,int>& NumPlaces, double time) {
        //std::cout << "[DEBUG] Applying gene regulation rules at time " << time << std::endl;
        for (size_t i = 0; i < geneRules.size(); ++i) {
            GeneRule &rule = geneRules[i];
            // Skip rule if already applied
            if (rule.applied) {
                //std::cout << "[DEBUG] GeneRule[" << i << "] already applied, skipping." << std::endl;
                continue;
            }
           // std::cout << "[DEBUG] Evaluating GeneRule[" << i << "]: Reaction = " << rule.reactionID;
            if (rule.timeSpecified) {
                //std::cout << ", Time condition: current time (" << time << ") >= " << rule.time;
            } else {
                //std::cout << ", Time condition: NA";
            }
            if (rule.placeSpecified && rule.thresholdSpecified) {
                //std::cout << ", Place condition: " << rule.place;
                if (NumPlaces.find(rule.place) != NumPlaces.end()) {
                    //int idx = NumPlaces[rule.place];
                    //double conc = Value[idx];
                    //std::cout << " (current conc = " << conc << ")";
                } else {
                    //std::cout << " (metabolite not found)";
                }
            } else {
               // std::cout << ", Place condition: NA";
            }
            //std::cout << ", compType = " << rule.compType << std::endl;

            bool triggered = true;
            if (rule.timeSpecified) {
                triggered = triggered && (time >= rule.time);
            }
            if (rule.placeSpecified && rule.thresholdSpecified) {
                if (NumPlaces.find(rule.place) != NumPlaces.end()) {
                    int idx = NumPlaces[rule.place];
                    double conc = Value[idx];
                    if (rule.compType == "min") {
                        triggered = triggered && (conc <= rule.threshold);
                    } else if (rule.compType == "max") {
                        triggered = triggered && (conc >= rule.threshold);
                    } else {
                        // Default behavior: use 'max' condition
                        triggered = triggered && (conc >= rule.threshold);
                    }
                } else {
                    triggered = false;
                }
            }
            //std::cout << "[DEBUG] GeneRule[" << i << "] triggered: " << (triggered ? "YES" : "NO") << std::endl;
            if (triggered) {
                for (size_t j = 0; j < vec_fluxb.size(); ++j) {
                    int colIdx = vec_fluxb[j].fromNametoid(rule.reactionID);
                    if (colIdx != -1) {
                        if (rule.newLB == rule.newUB) {
                            vec_fluxb[j].update_bound(colIdx, "GLP_FX", rule.newLB, rule.newUB);
                        } else {
                            vec_fluxb[j].update_bound(colIdx, "GLP_DB", rule.newLB, rule.newUB);
                        }
                        std::cout << "[DEBUG] Applied gene rule for reaction " << rule.reactionID 
                                  << " in LP problem " << j << " at time " << time 
                                  << " (new bounds: [" << rule.newLB << ", " << rule.newUB << "])" << std::endl;
                    } else {
                        std::cout << "[DEBUG] Reaction " << rule.reactionID 
                                  << " not found in LP problem " << j << std::endl;
                    }
                }
                rule.applied = true;
            }
        }
    }


    /**
     * Destructor to clean up memory allocated for the Vars pointer of pointers, ensuring no memory leaks.
     */
    ~FBAProcessor() {
        // Clean up memory for the array of pointers
        delete[] Vars;  // Only deallocates the array of pointers
        Vars = nullptr;  // Sets the pointer to nullptr to prevent invalid memory access
        FBGLPK::freeGLPKEnvironment();  // Sets the pointer to nullptr to prevent invalid memory access
        cout << "FBAProcessor destroyed and GLPK environment freed." << endl;
    }


		// -------------------------------------------------------------
		// λm  = 1 / Σk ( xB,k · xN,k · 10-12 )
		//       = 1e12 / Σk ( xB,k · xN,k )   [ xB in pgDW/cell ]
		//
		// moltiplicativeConstant  ≔ λm               (usato in updateFluxBounds)
		// -------------------------------------------------------------
		void calculateMultiplicativeConstant(double*                      Value,
				                                 const std::map<std::string,int>& NumPlaces,
				                                 std::size_t                   problemIndex,
				                                 const std::string&            reaction,
				                                 vector<class FBGLPK::LPprob>& vec_fluxb )
		{
				double denom_pgDWcell = 0.0;      // Σ xB,k · xN,k   (pgDW)
				const auto& lpSet =
				    reactionToFileMap.count(reaction) ?  reactionToFileMap.at(reaction)
				                                      : std::unordered_set<std::size_t>{problemIndex};

				std::cout << "\n[λ-DEBUG] reaction = " << reaction << '\n';
				std::cout << "model\tcells(xN)\txB(pgDW)\txB·xN\n";

				/* -------- ① somma di tutte le biomasse coinvolte -------- */
				for (auto prob : lpSet)
				{
				    long   xN = std::lround( std::floor(
				                 Value[ NumPlaces.at( problemBacteriaPlace.at(prob) ) ] ));
				    double xB = Value[ NumPlaces.at( problemBiomassPlace.at(prob) ) ];     // pgDW/cell

				    denom_pgDWcell += xB * static_cast<double>(xN);

				    std::cout << vec_fluxb[prob].getFilename() << '\t'
				              << xN   << '\t'
				              << xB   << '\t'
				              << xB * static_cast<double>(xN) << '\n';
				}

				/* -------- ② λ  -------- */
				if (denom_pgDWcell > 0.0)
				    multiplicativeConstant = trunc( 1e12 / denom_pgDWcell, decimalTrunc );
				else
				    multiplicativeConstant = 0.0;

				std::cout << "Σ xB·xN = " << denom_pgDWcell
				          << "  ⇒  λ = "  << multiplicativeConstant
				          << "  [mL/(gDW·h)]\n";
		}



    /**
     * Trims whitespace from both ends of a given string.
     * @param str The string to trim.
     * @return A trimmed string with no leading or trailing whitespace.
     */
    string trim(const string& str) {
        size_t first = str.find_first_not_of(' ');
        if (first == string::npos) return ""; // Return empty if string is made of spaces
        size_t last = str.find_last_not_of(' ');
        return str.substr(first, last - first + 1);
    }



    /**
     * Splits a string by a given delimiter and trims each resulting token.
     * @param str The string to split.
     * @param delim The delimiter character.
     * @return A vector of trimmed strings.
     */
    vector<string> splitAndTrim(const string& str, char delimiter) {
        vector<string> tokens;
        stringstream ss(str);
        string item;
        while (getline(ss, item, delimiter)) {
            item.erase(remove(item.begin(), item.end(), ' '), item.end()); // Remove spaces
            if (!item.empty()) {
                tokens.push_back(item);
            }
        }
        return tokens;
    }


// ── helper: to‑lower utility
static std::string toLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

// ── helper: build abbreviation (“cbd” from “Clostridium_butyricum_DSM_10702”)
static std::string modelAbbrev(const std::string& model)
{
    std::string abbr;
    bool newToken = true;
    for (char c : model) {
        if (std::isalpha(c)) {
            if (newToken) abbr += std::tolower(c);
            newToken = false;
        } else {
            newToken = true;          // hit '_' or digit etc. → next char starts a token
        }
    }
    return abbr;
}

// ──────────────────────────────────────────────────────────────
// smarter LP lookup: exact filename  ►  base‑name match  ►  abbreviation
// ──────────────────────────────────────────────────────────────
size_t findLPIndex(const std::vector<FBGLPK::LPprob>& vec_fluxb,
                                 const std::string& modelName) const
{
    // 1) exact match (unchanged behaviour)
    for (size_t i = 0; i < vec_fluxb.size(); ++i)
        if (vec_fluxb[i].getFilename() == modelName)
            return i;

    // 2) compare against base name (filename without extension, case‑insensitive)
    std::string target = toLower(modelName);
    for (size_t i = 0; i < vec_fluxb.size(); ++i) {
        std::string base = vec_fluxb[i].getFilename();
        base = base.substr(0, base.find_last_of('.'));   // strip extension
        if (toLower(base) == target) return i;
    }

    // 3) compare abbreviation with beginning of base name
    std::string abbr = modelAbbrev(modelName);
    for (size_t i = 0; i < vec_fluxb.size(); ++i) {
        std::string base = vec_fluxb[i].getFilename();
        base = base.substr(0, base.find_last_of('.'));
        base = toLower(base);
        if (base.find(abbr) == 0)        // base starts with abbreviation
            return i;
    }

    return std::numeric_limits<size_t>::max();           // not found
}


/**
 * @brief Ritorna il nome reazione finale, con logica speciale per Biomass.
 *
 * Se isBiomass == true:
 *   - se ha sia input che output => '_f'
 *   - se ha solo input => '_r'
 * Altrimenti, usa la logica di splitted e irreversibile.
 */
static std::string standardizeReactionName(
    const std::string& transitionName,
    const std::string& reactionBase,
    const std::vector<std::string>& inputs,
    const std::vector<std::string>& outputs,
    const std::map<std::string, unsigned int>& forwardMap,
    const std::map<std::string, unsigned int>& reverseMap,
    const std::map<std::string, unsigned int>& irreversMap,
    bool isBiomass
)
{
    std::cerr << "[DEBUG standardizeReactionName] "
              << "transition='"   << transitionName   << "', "
              << "reactionBase='"  << reactionBase     << "', "
              << "isBiomass="      << (isBiomass ? "true" : "false")
              << "\n";

    // 1) irreversibile puro
    if (irreversMap.find(reactionBase) != irreversMap.end()) {
        std::cerr << "  [DEBUG] '" << reactionBase
                  << "' è irreversibile => restituisco invariato\n";
        return reactionBase;
    }

    // 2) split disponibili?
    bool splittedF = (forwardMap.find(reactionBase + "_f") != forwardMap.end());
    bool splittedR = (reverseMap.find(reactionBase + "_r") != reverseMap.end());
    std::cerr << "  [DEBUG] splittedF=" << splittedF
              << ", splittedR=" << splittedR << "\n";

    // 3) nessun split => invariato
    if (!splittedF && !splittedR) {
        std::cerr << "  [WARN] '" << reactionBase
                  << "' non in irreversMap, forwardMap, né reverseMap\n";
        return reactionBase;
    }

    // 4) logica Biomass
    if (isBiomass) {
        std::cerr << "  [DEBUG Biomass] inputs=" << inputs.size()
                  << ", outputs=" << outputs.size() << "\n";
        if (!inputs.empty() && !outputs.empty() && splittedF) {
            std::cerr << "  [DEBUG Biomass] caso 4a → forward ('_f')\n";
            return reactionBase + "_f";
        }
        else if (!inputs.empty() && outputs.empty() && splittedR) {
            std::cerr << "  [DEBUG Biomass] caso 4b → reverse ('_r')\n";
            return reactionBase + "_r";
        }
        else {
            std::cerr << "  [WARN Biomass] mismatch, restituisco invariato\n";
            return reactionBase;
        }
    }

    // 5) logica “classica” non-biomass
    if (!outputs.empty() && splittedF) {
        std::cerr << "  [DEBUG] non-biomass con outputs → forward ('_f')\n";
        return reactionBase + "_f";
    }
    else if (!inputs.empty() && splittedR) {
        std::cerr << "  [DEBUG] non-biomass con inputs → reverse ('_r')\n";
        return reactionBase + "_r";
    }

    // 6) fallback
    std::cerr << "  [WARN] fallback finale → restituisco invariato\n";
    return reactionBase;
}


    /**
     * Reads and processes the FBA information from a specified file path, initializing data structures necessary for simulation.
     * This includes mapping each transition to its corresponding reaction, and associating reactions to LP problem indices.
     * Also builds a mapping from LP problem indices to sets of reactions for optimized problem solving.
     * @param filePath The path to the file containing FBA information.
     * @param vec_fluxb A vector of LP problems to be populated based on the file.
     * @return true if the file was read successfully, false otherwise.
     */
    bool readFBAInfo(const string& filePath, vector<class FBGLPK::LPprob>& vec_fluxb) {
        ifstream file(filePath);
        if (!file.is_open()) {
            cerr << "Failed to open file: " << filePath << endl;
            return false;
        }
        stringstream buffer;
        buffer << file.rdbuf();
        string content = buffer.str();
        file.close();

        size_t pos = 0;
        while ((pos = content.find('{', pos)) != string::npos) {
            size_t endPos = content.find('}', pos);
            if (endPos == string::npos) break;

            string jsonObject = content.substr(pos + 1, endPos - pos - 1);
            pos = endPos + 1;

            string transition, lpFile, reaction, multiplicityStr, bacteriaCountPlace, bacteriaBiomassPlace, inputPlacesStr, outputPlacesStr, isBiomassStr;
            int multiplicity = 1;
            bool isBiomass = false;

            size_t keyPos = 0, valuePos = 0;
            while ((keyPos = jsonObject.find('"', keyPos)) != string::npos) {
                size_t keyEnd = jsonObject.find('"', keyPos + 1);
                string key = jsonObject.substr(keyPos + 1, keyEnd - keyPos - 1);

                valuePos = jsonObject.find(':', keyEnd) + 1;
                size_t valueEnd = (jsonObject[valuePos] == '"') ? jsonObject.find('"', valuePos + 1) : jsonObject.find(',', valuePos);
                if (valueEnd == string::npos) valueEnd = jsonObject.length();
                string value = jsonObject.substr(valuePos, valueEnd - valuePos);
                value.erase(remove(value.begin(), value.end(), '"'), value.end());
                value.erase(remove(value.begin(), value.end(), ' '), value.end());

                if (key == "Transition") transition = value;
                else if (key == "LPfile") lpFile = value;
                else if (key == "Reaction") reaction = value;
                else if (key == "Multiplicity") multiplicityStr = value;
                else if (key == "BacteriaCountPlace") bacteriaCountPlace = value;
                else if (key == "BacteriaBiomassPlace") bacteriaBiomassPlace = value;
                else if (key == "InputPlaces") inputPlacesStr = value;
                else if (key == "OutputPlaces") outputPlacesStr = value;
                else if (key == "IsBiomass") isBiomass = (value == "true\n");

                keyPos = valueEnd + 1;
            }

            try {
                multiplicity = stoi(multiplicityStr);
            } catch (const std::invalid_argument& e) {
                cerr << "Error converting multiplicity to integer: " << multiplicityStr << endl;
                cerr << "Exception: " << e.what() << endl;
                return false;
            } catch (const std::out_of_range& e) {
                cerr << "Multiplicity value out of range: " << multiplicityStr << endl;
                cerr << "Exception: " << e.what() << endl;
                return false;
            }

            // Remove square brackets and split the strings
            inputPlacesStr.erase(remove(inputPlacesStr.begin(), inputPlacesStr.end(), '['), inputPlacesStr.end());
            inputPlacesStr.erase(remove(inputPlacesStr.begin(), inputPlacesStr.end(), ']'), inputPlacesStr.end());
            outputPlacesStr.erase(remove(outputPlacesStr.begin(), outputPlacesStr.end(), '['), outputPlacesStr.end());
            outputPlacesStr.erase(remove(outputPlacesStr.begin(), outputPlacesStr.end(), ']'), outputPlacesStr.end());

            auto removeUnwantedChars = [](string& str) {
                str.erase(remove_if(str.begin(), str.end(), [](char c) { return isspace(c) || iscntrl(c); }), str.end());
            };

            removeUnwantedChars(inputPlacesStr);
            removeUnwantedChars(outputPlacesStr);

            vector<string> inputs;
            if (!inputPlacesStr.empty() && !(inputPlacesStr == " ") && !(inputPlacesStr == "")) inputs = splitAndTrim(inputPlacesStr, ',');
            vector<string> outputs;
            if (!outputPlacesStr.empty() && !(outputPlacesStr == " ") && !(outputPlacesStr == "")) outputs = splitAndTrim(outputPlacesStr, ',');

            if (isBiomass) {
                if (!inputs.empty() && !outputs.empty()) {
                    biomassTransitions.insert(transition);
                    //std::cerr << "[DEBUG] => Inserito '" << transition << "' in biomassTransitions (entrambe input e output presenti).\n";
                } else {
                    //std::cerr << "[DEBUG] => '" << transition << "' è IsBiomass ma NON ha entrambi input e output => NON inserita in biomassTransitions.\n";
                }
            }

// … inside your while(jsonObject) loop, after you parse inputs, outputs, isBiomass, etc. …

string finalBaseRxn;
if (!lpFile.empty() && !reaction.empty()) {
    // 1) trova l'indice del modello LP
    size_t idx = findLPIndex(vec_fluxb, lpFile);
    if (idx == std::numeric_limits<size_t>::max()) {
        // file non trovato → nessuna FBA per questa transition
        transitionsWithoutFile.insert(transition);
    }
    else {
        // 2) normalizza con suffisso _f / _r
        finalBaseRxn = standardizeReactionName(
            transition,
            reaction,
            inputs,
            outputs,
            vec_fluxb[idx].getForwardReactions(),
            vec_fluxb[idx].getReverseReactions(),
            vec_fluxb[idx].getIrreversibileReactions(),
            isBiomass
        );

        // 3) verifica esistenza nella LP
        int col = vec_fluxb[idx].fromNametoid(finalBaseRxn);
        if (col < 0) {
            std::cerr << "[WARNING readFBAInfo] '" << finalBaseRxn
                      << "' non trovato in LP[" << vec_fluxb[idx].getFilename() << "]\n";
            return false;
        }

        // 4) registra metadata GLPK‑side
        problems.insert(idx);
        FBAproblems[transition] = idx;
        problemsToReactions[idx].insert(transition);
        bacteriaToBioMax[idx]  = vec_fluxb[idx].getBioMax();
        bacteriaToBioMean[idx] = vec_fluxb[idx].getBioMean();
        bacteriaToBioMin[idx]  = vec_fluxb[idx].getBioMin();

        if (!bacteriaCountPlace.empty() && bacteriaCountPlace != "N/A") {
            problemBacteriaPlace[idx]       = bacteriaCountPlace;
            reactionToBacteria[transition]  = bacteriaCountPlace;
            hasMultiSpecies = true;
        }
        if (!bacteriaBiomassPlace.empty() && bacteriaBiomassPlace != "N/A") {
            problemBiomassPlace[idx]           = bacteriaBiomassPlace;
            reactionToBacteriaBIOMASS[transition] = bacteriaBiomassPlace;
            hasBioMASS = true;
        }

        // 5) genera gli ID ESPN (uno per specie se biomass, altrimenti uguale al base)
        if (isBiomass) {
            // estrai suffisso specie da "…_in_cbd1" → "cbd1"
            string species = transition.substr(transition.find_last_of('_') + 1);
            // UNIQUE ESPN reaction ID: include _f/_r e specie
            string uniqRxn = finalBaseRxn + "_" + species;

            // ESPN‑side mapping
            FBAreact[transition]    = uniqRxn;
            FBAreactions.insert(uniqRxn);
            // il posto corrispondente è la biomassa di quella specie
            FBAplace[uniqRxn]       = problemBiomassPlace[idx];
            // salva la corrispondenza per tornare al base‑rxn GLPK
            speciesToBaseRxn[uniqRxn] = finalBaseRxn;

            std::cout << "[DEBUG readFBAInfo] transition=" << transition
                      << "  uniqRxn=" << uniqRxn
                      << "  FBAplace[" << uniqRxn << "]=" << FBAplace[uniqRxn]
                      << "\n";
        }
        else {
            // non‑biomass: ESPN ID = baseReaction (_f/_r incluso)
            FBAreact[transition]    = finalBaseRxn;
            FBAreactions.insert(finalBaseRxn);

            // merge inputs+outputs senza duplicati
            vector<string> newPlaces = inputs;
            newPlaces.insert(newPlaces.end(), outputs.begin(), outputs.end());
            vector<string> allPlaces = splitAndTrim(FBAplace[finalBaseRxn], ',');
            for (auto &p : newPlaces) {
                if (find(allPlaces.begin(), allPlaces.end(), p) == allPlaces.end()) {
                    allPlaces.push_back(p);
                }
            }
            // ricomponi lista
            string merged;
            for (size_t i = 0; i < allPlaces.size(); ++i) {
                if (i) merged += ",";
                merged += allPlaces[i];
            }
            FBAplace[finalBaseRxn] = merged;

            std::cout << "[DEBUG readFBAInfo] transition=" << transition
                      << "  finalBaseRxn=" << finalBaseRxn
                      << "  FBAplace[" << finalBaseRxn << "]={" << merged << "}\n";
        }
    }
}


// finally, multiplicity and gene rules wiring
ReactionMultiplicity[transition] = multiplicity;
if (!bacteriaCountPlace.empty() && bacteriaCountPlace != "N/A") {
    reactionToBacteria[transition]        = bacteriaCountPlace;
}
if (!bacteriaBiomassPlace.empty() && bacteriaBiomassPlace != "N/A") {
    reactionToBacteriaBIOMASS[transition] = bacteriaBiomassPlace;
}


        }

        return true;
    }


		double dynamicThreshold(double previousConcentration) {
    		double baseThreshold = 1e-3;
    		return std::max(minDeltaThreshold, baseThreshold * fabs(previousConcentration));
    		// return minDeltaThreshold;
		}

		/**
		 * Verifica quali metaboliti hanno subito un cambiamento significativo
		 * e, per ciascuno, stampa un log dettagliato delle condizioni valutate.
		 */
		set<string> checkSignificantChange(double* Value,
				                               map<string, int>& NumPlaces,
				                               double time)
		{
				// Pulisco la lista delle reazioni da aggiornare
				reactionsToUpdate.clear();

				// Il set che restituiremo
				set<string> changedMetabolites;

				// Ricavo la soglia relativa effettiva
				double relChange = (relEpsilon > 100.0) 
				                   ? 100.0 
				                   : ((relEpsilon < 0.0 && absEpsilon == -1) 
				                       ? 0 
				                       : ((relEpsilon < 0.0 && absEpsilon != -1) 
				                           ? -1 
				                           : relEpsilon));

				bool useAbsolute = (absEpsilon != -1);
				bool useRelative = (relChange != -1);

				// Ciclo su tutti i metaboliti coinvolti nelle FBA
				for (const auto& place : FBAplace) {
				    // DEBUG: mostro la reaction e la stringa dei posti
				    std::cout << "[DBG check] reaction=" << place.first
				              << "  placesStr=\"" << place.second << "\"\n";

				    // suddivido la stringa per ottenere singoli nomi di metaboliti
				    auto mets = splitAndTrim(place.second, ',');
				    for (auto& metabolite : mets) {
				        // DEBUG: nome del metabolita corrente
				        std::cout << "[DBG check]   metabolite=\"" << metabolite << "\"\n";

				        // controllo che esista in NumPlaces
				        auto itNP = NumPlaces.find(metabolite);
				        if (itNP == NumPlaces.end()) {
				            std::cerr << "[ERROR check]   NumPlaces ha NO entry per \"" 
				                      << metabolite << "\" → skip\n";
				            continue;
				        }

				        // Stato corrente e precedente (troncamento incluso)
				        double currentConcentration  = trunc(Value[itNP->second], decimalTrunc);
				        double previousConcentration = previousConcentrations[metabolite];

				        // Calcolo differenze assoluta e percentuale
				        double absoluteDiff  = fabs(currentConcentration - previousConcentration);
				        double percentChange = (previousConcentration != 0.0)
				                               ? (absoluteDiff / fabs(previousConcentration)) * 100.0
				                               : 100.0;

				        // Valuto le quattro condizioni
				        bool condAbsolute = useAbsolute && (absoluteDiff  > absEpsilon);
				        bool condRelative = useRelative && (percentChange > relChange);
				        bool condFromZero = useRelative && (previousConcentration == 0.0 && currentConcentration != 0.0);
				        bool condFirstSol = (previousConcentration == -1.0);

				        // Debug: stampo tutti i valori e l’esito di ciascuna condizione
				        std::cout
				            << "[DBG] met="   << metabolite
				            << " prev="       << previousConcentration
				            << " curr="       << currentConcentration
				            << " Δabs="       << absoluteDiff
				            << " %Δ="         << percentChange << "%"
				            << " -> abs? "    << (condAbsolute  ? "YES" : "no")
				            << " rel? "       << (condRelative  ? "YES" : "no")
				            << " zero? "      << (condFromZero  ? "YES" : "no")
				            << " first? "     << (condFirstSol  ? "YES" : "no")
				            << std::endl;

				        // Se almeno una condizione è vera, segno il metabolita per l’aggiornamento
				        if (condAbsolute || condRelative || condFromZero || condFirstSol) {
				            std::cout << "      >> will update metabolite: " << metabolite << std::endl;
				            changedMetabolites.insert(metabolite);
				        }
				    }
				}

				return changedMetabolites;
		}



    /**
     * Truncates a double value to a specified number of decimal places.
     * This method uses mathematical manipulation to achieve truncation without rounding.
     * The value is scaled up by a power of ten corresponding to the desired number of decimal places,
     * then truncated using floor to drop any fractional part, and finally scaled back down.
     * This method ensures that the truncation is always downward, similar to how floor operates on positive numbers.
     *
     * @param value The double value to truncate.
     * @param decimal The number of decimal places to retain in the truncated result.
     * @return double The truncated value, preserving only the specified number of decimals.
     * 
     * Example:
     *    If value = 123.456789 and decimal = 3,
     *    the result would be 123.456.
     */
    double trunc(double value, double decimal) {
        const double multiplier = pow(10.0, decimal); // Calculate the scaling factor based on the number of decimals
        return floor(value * multiplier) / multiplier; // Scale up, truncate, scale down
        // Alternative implementation using type casting:
        // return ((unsigned long int)(value * multiplier)) / multiplier; // Cast to integer to truncate
    }


	void updateFluxBounds(
    double* Value,
    std::vector<FBGLPK::LPprob>& vec_fluxb,
    std::map<std::string,int>& NumPlaces,
    std::set<std::string>& changedMetabolites,
    double time
) {
    std::set<size_t> problemsToUpdate;

    std::cout << "[DEBUG updateFluxBounds] time=" << time 
              << ", changedMetabolites={";
    for (const auto &m : changedMetabolites) std::cout << m << ",";
    std::cout << "}\n";

    if (changedMetabolites.empty()) {
        // caso: cambiata solo biomassa → aggiorno comunque tutte le EX_
        std::cout << "[DEBUG updateFluxBounds] No metabolites → forcing EX_ updates\n";
        for (const auto& rxn : reactionsToUpdate) {
            auto itMap = reactionToFileMap.find(rxn);
            if (itMap == reactionToFileMap.end()) continue;
            for (auto prob : itMap->second) {
                bool alive = (deadBacterialSpecies.find(prob) == deadBacterialSpecies.end()
                              || !deadBacterialSpecies[prob]);
                std::cout << "    LP#" << prob 
                          << (alive ? " alive" : " dead") << "\n";
                if (alive) problemsToUpdate.insert(prob);
            }
        }
    }
    else {
        // caso normale: solo problemi legati ai metaboliti cambiati
        for (const auto& metabolite : changedMetabolites) {
            std::cout << "[DEBUG updateFluxBounds] checking metabolite → " << metabolite << "\n";
            auto itProb = metaboliteToProblems.find(metabolite);
            if (itProb == metaboliteToProblems.end()) {
                std::cout << "    no LP problems mapped for " << metabolite << "\n";
                continue;
            }
            for (auto prob : itProb->second) {
                bool alive = (deadBacterialSpecies.find(prob) == deadBacterialSpecies.end()
                              || !deadBacterialSpecies[prob]);
                std::cout << "    LP#" << prob 
                          << (alive ? " alive" : " dead") << "\n";
                if (alive) problemsToUpdate.insert(prob);
            }
        }
    }

    // 2) per ciascun problema, aggiorno le reazioni in reactionsToUpdate
    for (auto problemIndex : problemsToUpdate) {
        std::cout << "[DEBUG updateFluxBounds] Solving LP #" << problemIndex
                  << " (" << vec_fluxb[problemIndex].getFilename() << ")\n";

        for (const auto& reaction : reactionsToUpdate) {
            std::cout << "  [DEBUG] reaction=" << reaction << "\n";

            // skip biomass reactions (tutte quelle che iniziano con "EX_bio")
            if (reaction.rfind("EX_bio", 0) == 0) {
                std::cout << "    skipped: biomass transition\n";
                continue;
            }

            int colIdx = vec_fluxb[problemIndex].fromNametoid(reaction);
            if (colIdx < 0) {
                std::cout << "    skipped: not in LP\n";
                continue;
            }
            // solo reverse (_r)
            if (!endsWith(reaction, "_r")) {
                std::cout << "    skipped: not a reverse (_r)\n";
                continue;
            }

            // estraggo la stringa dei posti da FBAplace
            auto itFP = FBAplace.find(reaction);
            if (itFP == FBAplace.end()) {
                std::cerr << "    [ERROR] no FBAplace entry for reaction '" 
                          << reaction << "'\n";
                continue;
            }
            const std::string& placesStr = itFP->second;
            std::cout << "    [DEBUG] FBAplace['" << reaction 
                      << "'] = \"" << placesStr << "\"\n";

            // split e somma le concentrazioni
            auto mets = splitAndTrim(placesStr, ',');
            if (mets.empty()) {
                std::cerr << "    [ERROR] no metabolites parsed from placesStr\n";
                continue;
            }

            double pooledConc = 0.0;
            for (auto& met : mets) {
                std::cout << "    [DEBUG] handling metabolite place → " << met << "\n";
                auto itNP = NumPlaces.find(met);
                if (itNP == NumPlaces.end()) {
                    std::cerr << "      [ERROR] NumPlaces has no entry for '" 
                              << met << "'\n";
                    continue;
                }
                double conc = trunc(Value[itNP->second], decimalTrunc);
                std::cout << "      C_m(" << met << ") = " << conc << "\n";
                pooledConc += conc;
            }
            std::cout << "    [DEBUG] pooled C_m = " << pooledConc << "\n";

            // calcolo λ e vr_max_h
            calculateMultiplicativeConstant(Value, NumPlaces, problemIndex, reaction, vec_fluxb);
            double lambda    = multiplicativeConstant;
            double vr_max_h  = systemExternalVolume * pooledConc * lambda;
            std::cout << "    λ = " << lambda 
                      << ", vr_max_h = " << vr_max_h << "\n";

            if (pooledConc <= 0.0 || lambda == 0.0) {
                vr_max_h = 0.0;
                std::cout << "    adjusted vr_max_h to zero due to conc or λ\n";
            }

            // applico il bound
            double lb = 0.0;
            const char* type = (vr_max_h <= 0.0 ? "GLP_FX" : "GLP_DB");
            vec_fluxb[problemIndex].update_bound(
                colIdx,
                type,
                lb,
                trunc(vr_max_h, decimalTrunc)
            );
            std::cout << "    update_bound(col=" << colIdx
                      << ", type=" << type
                      << ", lb=" << lb
                      << ", ub=" << trunc(vr_max_h, decimalTrunc)
                      << ")\n";
        }
    }
}



    /**
     * Solves the FBA linear programming problems that are affected by significant changes in metabolite concentrations.
     * This method iterates through a subset of LP problems identified by changed metabolites, solving them to find the optimal flux
     * distributions under current constraints. It updates an array of variable pointers, storing the solution
     * variables for each problem. This is crucial for subsequent steps in metabolic analysis where these
     * solutions are used to compute metabolic flux rates.
     * 
     * This method improves efficiency by focusing on problems affected by recent changes, ensuring that
     * results are immediately accessible and the model is adjusted only as needed in response to
     * changing conditions within the biological system.
     * 
     * @param vec_fluxb Reference to a vector containing FBGLPK::LPprob objects.
     * @param changedMetabolites Set of metabolites that have changed, used to identify relevant LP problems to solve.
     */
    void solveFBAProblems(vector<class FBGLPK::LPprob>& vec_fluxb, set<string>& changedMetabolites) {
        set<size_t> toSolve;
        for (const string& metabolite : changedMetabolites) {
            if (metaboliteToProblems.find(metabolite) != metaboliteToProblems.end()) {
                for (auto problemIndex : metaboliteToProblems[metabolite]) {
                    if ((deadBacterialSpecies.find(problemIndex) == deadBacterialSpecies.end() || !deadBacterialSpecies[problemIndex])) {
                        toSolve.insert(problemIndex);
                    }
                }
            }
        }

        for (auto index : toSolve) {
            //cout << "Problema da risolvere ---> " << index <<  "   Name: " << vec_fluxb[index].getFilename() << endl;
            vec_fluxb[index].solve();
            // PARSIMONIUS FLAGS:
            if(vec_fluxb[index].getPFbaFlag() != -1) performPFBA(vec_fluxb, index);
            Vars[index] = vec_fluxb[index].getVariables(); // Update the variable pointers
        }
    }


    /**
     * Updates the stored concentrations of metabolites based on the latest computational results.
     * This method iterates over a map linking metabolite names to their indices in an array of current
     * metabolic concentrations. It updates a map of previous concentrations with the current values,
     * which is essential for tracking changes in metabolite levels over time and responding to dynamic
     * metabolic conditions. This tracking supports the system's ability to determine significant changes
     * in metabolic states that may require further adjustments in the model.
     * 
     * By maintaining an up-to-date record of metabolite concentrations, this method ensures that
     * the metabolic simulation reflects the most current state of the system, enabling accurate
     * and timely decision-making in metabolic engineering and research.
     * 
     * @param Value Pointer to an array containing the current concentrations of metabolites, indexed numerically.
     * @param NumPlaces Map linking metabolite names to their respective indices in the Value array, used for direct access.
     * @param changedMetabolites Set of metabolite names whose concentrations have changed significantly.
     */
    void updateConcentrations(double* Value, map<string, int>& NumPlaces, set<string>& changedMetabolites) {
        for (const string& metabolite : changedMetabolites) {
            int index = NumPlaces.at(metabolite);
            double newConcentration = trunc(Value[index], decimalTrunc);
            previousConcentrations[metabolite] = newConcentration; // Update the concentration record even if negative
            std::cout << "[DBG updateConcentrations] " 
						    << metabolite 
						    << " ← " << newConcentration 
						    << std::endl;

        }
    }


		/**
		 * Updates the flux bounds and solves the FBA problems only per le reazioni 
		 * realmente collegate ai metaboliti che hanno subito un cambiamento significativo.
		 *
		 * @param Value Pointer to l’array di valori metabolici correnti.
		 * @param vec_fluxb Vector di oggetti LPprob.
		 * @param NumPlaces Mappa place→indice.
		 * @param time     Tempo corrente (passato per compatibilità, ma non usato qui).
		 */
			void updateFluxBoundsAndSolve(
					double* Value,
					std::vector<FBGLPK::LPprob>& vec_fluxb,
					std::map<std::string,int>& NumPlaces,
					double time
			) {
					// 1) individua tutti i metaboliti cambiati
					std::set<std::string> origChanged = checkSignificantChange(Value, NumPlaces, time);
					if (origChanged.empty()) {
						  std::cout << "[DEBUG updateFluxBoundsAndSolve] No metabolites changed (incl. biomassa) → still updating FBAplace reactions\n";
					}

					// 2) escludi i posti di biomassa
					std::set<std::string> changed = origChanged;
					for (const auto& kv : problemBiomassPlace) {
						  const std::string& bioPlace = kv.second;
						  if (changed.erase(bioPlace)) {
						      std::cout << "[DEBUG updateFluxBoundsAndSolve] Excluding biomass place '"
						                << bioPlace << "' from dynamic update\n";
						  }
					}

					// 3) se è cambiata solo la biomassa, forziamo comunque l’aggiornamento
					if (!origChanged.empty() && changed.empty()) {
						  std::cout << "[DEBUG updateFluxBoundsAndSolve] Only biomass places changed → forcing FBAplace reactions update\n";
					}

					// 4) costruisco il set di reazioni da aggiornare:
					reactionsToUpdate.clear();
					//   a) quelle legate ai metaboliti effettivamente cambiati
					for (const auto& met : changed) {
						  auto it = metaboliteToReactions.find(met);
						  if (it != metaboliteToReactions.end()) {
						      for (const auto& rxn : it->second) {
						          reactionsToUpdate.insert(rxn);
						      }
						  }
					}
					//   b) tutte le reazioni che compaiono in FBAplace (quindi realmente FBA), 
					//      tranne quelle di biomassa (prefisso "EX_bio")
					for (const auto& kv : FBAplace) {
						  const std::string& rxn = kv.first;
						  if (rxn.rfind("EX_", 0) == 0 && rxn.rfind("EX_bio", 0) != 0) {
						      reactionsToUpdate.insert(rxn);
						  }
					}

					std::cout << "[DEBUG updateFluxBoundsAndSolve] Updating "
						        << reactionsToUpdate.size()
						        << " reactions for "
						        << (changed.empty() ? "0" : std::to_string(changed.size()))
						        << " metabolites\n";

					// 5) applica eventuali regole di gene regulation
					applyGeneRegulationRules(Value, vec_fluxb, NumPlaces, time);

					// 6) aggiorna i bound di flusso per tutte le reactionsToUpdate
					updateFluxBounds(Value, vec_fluxb, NumPlaces, changed, time);
					std::cout << "[DEBUG updateFluxBoundsAndSolve] Finished updateFluxBounds\n";

					// 7) aggiorna i bound non-FBA (_r e _f interni)
					updateNonFBAReactionsUB(vec_fluxb, NumPlaces, Value);
					std::cout << "[DEBUG updateFluxBoundsAndSolve] Finished updateNonFBAReactionsUB\n";

					// 8) risolvi i problemi FBA: se changed è vuoto, uso origChanged per il solve
					const std::set<std::string>& toSolveMets = changed.empty() ? origChanged : changed;
					solveFBAProblems(vec_fluxb, const_cast<std::set<std::string>&>(toSolveMets));
					std::cout << "[DEBUG updateFluxBoundsAndSolve] Finished solveFBAProblems\n";

					// 9) aggiorna le concentrazioni precedenti solo per i metaboliti non-biomassa cambiati
					updateConcentrations(Value, NumPlaces, changed);
					std::cout << "[DEBUG updateFluxBoundsAndSolve] Updated previous concentrations\n";

					count += 1;
			}



			/***************************************************************************
			 *  Compute reaction-rate (INTENSITY) for an ESPN transition                *
			 *  - segue le definizioni:                                                *
			 *      - Eq.(\ref{eq:zeta_biomass})  →  ζ_B = MW_B · x_B  (MW_B = 1)      *
			 *      - Eq.(\ref{eq:zeta_non_biomass}) → ζ(m) = x_N · x_B · 10⁻¹²        *
			 *  - v★  è l’output FBA (mmol / gDW / h)                                   *
			 *  - rate restituito:                                                     *
			 *      - biomassa  →  pgDW / cell / h  (può essere usato per Δx_B)        *
			 *      - altri metaboliti → mmol / h                                       *
			 ***************************************************************************/
			double computeRate(std::vector<FBGLPK::LPprob>&  vec_fluxb,
						             std::map<std::string,int>&    NumPlaces,
						             const std::vector<std::string>& NameTrans,
						             const double*                 Value,
						             int                           T,
						             double                        decimalTrunc,
						             const double&                 time)
			{
					const std::string& transitionName = NameTrans[T];
					const size_t       problemIndex   = FBAproblems[transitionName];
    // 1) ricavo l'ID ESPN (con specie+_f/_r)
					std::string espnRxn = FBAreact[transitionName];
					// 2) se è biomass, mappalo al nome base GLPK
					std::string glpkRxn = speciesToBaseRxn.count(espnRxn)
						                   ? speciesToBaseRxn[espnRxn]
						                   : espnRxn;

					std::cout << "\n[computeRate] transition=" << transitionName
						        << "  espnRxn=" << espnRxn
						        << "  glpkRxn=" << glpkRxn << std::endl;
					/*―― 1. ricavo il flusso specifico ottimo v★ [mmol/gDW/h] ――*/
					int varIndex   = vec_fluxb[problemIndex].fromNametoid(glpkRxn);
					double vStar   = Vars[problemIndex][varIndex];               // mmol / gDW / h
					double mult    = ReactionMultiplicity[transitionName];       // fattore d’unità (ad es. mM→mmol)

					/*―― 2. stato attuale sistema ――*/
					double xB = 1.0;   // pgDW / cell
					if (hasBioMASS)
						  xB = Value[ NumPlaces.at( problemBiomassPlace.at(problemIndex) ) ];   // x_B

					long   xN = 1;     // numero di cellule
					if (hasMultiSpecies)
						  xN = std::lround( std::floor( Value[ NumPlaces.at( problemBacteriaPlace.at(problemIndex) ) ] ) );

					const bool isBiomass = (biomassTransitions.find(transitionName) != biomassTransitions.end());

					/*―― 3. calcolo ζ(m) secondo il tipo di transizione ――*/
					double zeta = 1.0;
					double rawZeta = 1.0;
					if (isBiomass) {
						  /*  ζ_B  =  MW_B · x_B   (MW_B = 1 gDW/mmol)                       */
						  zeta = xB;                                    // pgDW / cell
					} else {
							// ζ(m) = (x_N · x_B · 10⁻¹²) / V_ext   [gDW per unit volume]
							rawZeta = static_cast<double>(xN) * xB * 1e-12;
							zeta = rawZeta / systemExternalVolume;
					}

					/*―― 4. rate finale ――*/
					double rate = vStar * zeta * mult;
									
					/*―― 5. DEBUG VERBOSO ――*/
					std::cout << "\n[computeRate] transition=" << transitionName
										<< (isBiomass ? "  (biomass)" : "  (boundary)")
										<< "\n   v★                  = " << vStar                        << "  [mmol/gDW/h]"
										<< "\n   x_B                 = " << xB                           << "  [pgDW/cell]"
										<< "\n   x_N                 = " << xN                           << "  [cell]"
										<< "\n   raw ζ (x_N·x_B·10⁻¹²)= " << rawZeta                      << "  [gDW]"
										<< "\n   V_ext               = " << systemExternalVolume         << "  [volume units]"
										<< "\n   ζ_scaled (raw ζ/V_ext)= " << zeta                        << "  [gDW per vol]"
										<< "\n   multiplicity (μ)     = " << mult                         << ""
										<< "\n   → rate              = " << rate
										<< (isBiomass ? "  [pgDW/cell/h]" : "  [mmol/(h·vol)]")
										<< std::endl;

					return rate;
			}

		// Updated initializeConcentrations and mapMetabolitesToProblems

		/**
		 * Initializes previousConcentrations for each unique metabolite
		 * (splitting comma-separated lists in FBAplace).
		 */
		void initializeConcentrations(double* Value, map<string, int>& NumPlaces) {
				// Collect all unique metabolites
				unordered_set<string> allMets;
				for (const auto& kv : FBAplace) {
				    auto mets = splitAndTrim(kv.second, ',');
				    for (const auto& met : mets) {
				        allMets.insert(met);
				    }
				}

				// Initialize previousConcentrations once per metabolite
				previousConcentrations.clear();
				cout << "[DEBUG initializeConcentrations] Initializing "
				     << allMets.size() << " metabolites" << endl;
				for (const auto& met : allMets) {
				    previousConcentrations[met] = -1.0;
				    cout << "  metabolite '" << met << "' ← -1" << endl;
				}
		}

		/**
		 * Maps each unique metabolite to the set of LP problems that it affects.
		 * Splits comma-separated lists in FBAplace and avoids duplicate mappings.
		 */
		void mapMetabolitesToProblems() {
				metaboliteToProblems.clear();
				cout << "[DEBUG mapMetabolitesToProblems] Mapping metabolites to LP problems" << endl;
				
				for (const auto& kv : FBAproblems) {
				    const string& transition = kv.first;
				    size_t        pidx       = kv.second;
				    const string& reaction   = FBAreact[transition];
				    auto it = FBAplace.find(reaction);
				    if (it == FBAplace.end()) continue;

				    auto mets = splitAndTrim(it->second, ',');
				    for (const auto& met : mets) {
				        auto& probSet = metaboliteToProblems[met];
				        if (probSet.insert(pidx).second) {
				            cout << "  metabolite '" << met
				                 << "' → LP problem #" << pidx << endl;
				        }
				    }
				}

				// Final summary of mappings
				cout << "[DEBUG mapMetabolitesToProblems] Final summary:" << endl;
				for (const auto& kv : metaboliteToProblems) {
				    cout << "  '" << kv.first << "' → { ";
				    for (auto p : kv.second) cout << p << " ";
				    cout << "}" << endl;
				}
		}


    /**
     * Initializes the data structures necessary for the FBA simulation process.
     * This method is invoked only once per instance, conforming to the singleton design pattern of this class.
     * It loads necessary data from a specified .fbainfo file into various structures for managing reactions,
     * places, and problem indices, and initializes all metabolic concentrations to a default negative value.
     * This ensures that any uninitialized metabolic state is clearly indicated and handled appropriately.
     *
     * @param NameTrans The names of the transitions involved in the network, used for initialization checks.
     * @param vec_fluxb A reference to a vector of LP problems representing the metabolic fluxes.
     * @param Value Pointer to an array storing the metabolic concentrations.
     * @param NumPlaces A map linking place names to their corresponding indices.
     * @return void
     * @throws If the file cannot be read, it prints an error message and terminates initialization.
     */
    void init_data_structures_class(const vector<string>& NameTrans, vector<class FBGLPK::LPprob>& vec_fluxb, double* Value, map<string, int>& NumPlaces) {
        if (!readFBAInfo(FBA_INFO_FILE, vec_fluxb)) { // Attempt to read initialization data from the .fbainfo file
            cerr << "Failed to read places from .fbainfo file" << endl;
            exit(1); // Early return if file reading fails
        }

        Vars = new double*[vec_fluxb.size()];
        initializeConcentrations(Value, NumPlaces);
        mapMetabolitesToProblems();
    }


		void mapReactionsFromProblems(const vector<class FBGLPK::LPprob>& vec_fluxb) {
    		reactionToFileMap.clear(); // Clear the existing map
    		for (const auto& index : problems) {
        		const auto& lp = vec_fluxb[index];
        		vector<string> reactions = lp.getReactions();
        		for (const auto& reaction : reactions) {
            		if (reaction.substr(0, 3) == "EX_") { // Check if reaction starts with 'EX_'
                		reactionToFileMap[reaction].insert(index);
            		}
        		}
    		}
		}

		// -------------------------------------------------------------------------
		// Utility: choose first existing filename from a list                   
		// -------------------------------------------------------------------------
		static std::string firstExisting(const std::vector<std::string>& paths)
		{
				for (const auto& p : paths) {
				    std::ifstream test(p.c_str());
				    if (test.good()) return p;
				}
				return "";            // none found
		}

		// Unified loader for non-projected (_f and _r) bounds from a single CSV with debug
		// For forward _f: apply the baseUb directly to the LP; for reverse _r: only store baseUb for runtime updates
		void updateNonFBAReactionUpperBoundsFromFile(
				std::vector<FBGLPK::LPprob>& vec_fluxb,
				const std::string& /*unused*/)
		{
				std::cout << "[DEBUG NonFBARead] Starting non-FBA bounds loading\n";

				// Trova il CSV
				const std::string path = firstExisting({
				    "non_projected_bounds_gui.csv",
				    "non_projected_bounds.csv",
				    "EX_upper_bounds_nonFBA.csv"
				});
				if (path.empty()) {
				    std::cerr << "[Warning] non-FBA bounds CSV not found – UBs unchanged.\n";
				    return;
				}
				std::cout << "[DEBUG NonFBARead] Using file: " << path << "\n";

				std::ifstream file(path);
				if (!file.good()) {
				    std::cerr << "[Warning] Cannot open non-FBA bounds file '" << path
				              << "' – UBs unchanged.\n";
				    return;
				}

				using Key = std::pair<std::string, std::size_t>;
				using DebugRow = std::tuple<std::string, std::string, double, double, double>;
				std::unordered_set<Key, RxnLPHash> done;
				std::vector<DebugRow> debugRows;

				// Salta header
				std::string header;
				std::getline(file, header);
				std::cout << "[DEBUG NonFBARead] Header: " << header << "\n";

				std::string line;
				while (std::getline(file, line)) {
				    std::cout << "[DEBUG NonFBARead] Raw line: " << line << "\n";
				    if (line.empty()) continue;
				    std::istringstream ss(line);
				    std::string reaction, model, bgStr, volStr;
				    std::getline(ss, reaction, ',');
				    std::getline(ss, model,    ',');
				    std::getline(ss, bgStr,    ',');
				    std::getline(ss, volStr,   ',');
				    std::cout << "[DEBUG NonFBARead] Parsed: reaction='" << reaction
				              << "', model='" << model
				              << "', background='" << bgStr
				              << "', volume='" << volStr << "'\n";
				    if (reaction.empty() || bgStr.empty() || volStr.empty()) {
				        std::cout << "[DEBUG NonFBARead] Skipping incomplete line\n";
				        continue;
				    }

				    double background = 0.0, volume = 1.0;
				    try {
				        background = std::stod(bgStr);
				        volume     = std::stod(volStr);
				    } catch (const std::exception& e) {
				        std::cerr << "[NonFBA] Warning: non-numeric background/volume for '"
				                  << reaction << "': " << e.what() << "\n";
				        continue;
				    }

				    double baseUb = fabs(background) * volume;
				    std::cout << "[DEBUG NonFBARead] baseUb = |" << background << "| * "
				              << volume << " = " << baseUb << "\n";

				    // Trova gli LP corrispondenti
				    std::vector<std::size_t> lpList;
				    std::size_t idx = findLPIndex(vec_fluxb, model);
				    if (idx != std::numeric_limits<std::size_t>::max()) {
				        lpList.push_back(idx);
				        std::cout << "[DEBUG NonFBARead] Model matched LP#" << idx
				                  << " (" << vec_fluxb[idx].getFilename() << ")\n";
				    }
				    else if (reactionToFileMap.count(reaction)) {
				        for (auto p : reactionToFileMap[reaction]) {
				            lpList.push_back(p);
				            std::cout << "[DEBUG NonFBARead] reactionToFileMap → LP#" << p
				                      << "\n";
				        }
				    }

				    if (lpList.empty()) {
				        std::cerr << "[NonFBA] Warning: no LP found for '" << reaction << "'.\n";
				        continue;
				    }

				    // Determina subtype: external se la parte prima di _r/_f contiene "_e"
				    bool isExternal = false;
				    if (endsWith(reaction, "_r") || endsWith(reaction, "_f")) {
				        std::string prefix = reaction.substr(0, reaction.size() - 2);
				        if (prefix.find("_e") != std::string::npos) {
				            isExternal = true;
				        }
				    }
				    std::string subtype = isExternal ? "exchange" : "internal";
				    std::cout << "[DEBUG NonFBARead] reaction '" << reaction
				              << "' → subtype = " << subtype << "\n";

				    // Registra baseUb e subtype per ogni LP, sia _f che _r
				    for (auto probIdx : lpList) {
				        Key key{reaction, probIdx};
				        if (done.count(key)) {
				            std::cout << "[DEBUG NonFBARead] Already registered key for '"
				                      << reaction << "' LP#" << probIdx << "\n";
				            continue;
				        }
				        done.insert(key);
				        NonFBAReactionBaseUB[key] = baseUb;
				        ReactionSubtype[key]      = subtype;
				        std::cout << "[DEBUG NonFBARead] Stored: key=(" << reaction
				                  << ", LP#" << probIdx << "), baseUb=" << baseUb
				                  << ", subtype=" << subtype << "\n";

				        debugRows.emplace_back(
				            reaction,
				            vec_fluxb[probIdx].getFilename(),
				            background,
				            volume,
				            baseUb
				        );
				    }
				}
				file.close();
				std::cout << "[DEBUG NonFBARead] Finished parsing CSV, dumping debug CSV\n";
/*
				// Genera il CSV di debug
				if (!debugRows.empty()) {
				    std::ofstream dbg("debug_nonFBA_bounds.csv");
				    dbg << "reaction,LPfile,background_conc,volume,baseUb\n";
				    for (auto &r : debugRows) {
				        dbg << std::get<0>(r) << ','
				            << std::get<1>(r) << ','
				            << std::get<2>(r) << ','
				            << std::get<3>(r) << ','
				            << std::get<4>(r) << '\n';
				    }
				    dbg.close();
				    std::cout << "[DEBUG NonFBARead] Debug CSV written: debug_nonFBA_bounds.csv\n";
				} else {
				    std::cout << "[DEBUG NonFBARead] No debug rows to write\n";
				}*/
		}
		
		/**
		 * Adapted loader for projected FBA bounds from a single CSV (handles only forward _f; skips reverse _r)
		 * Also captures systemExternalVolume from the first data row and seeds NonFBAReactionBaseUB & ReactionSubtype
		 * so that updateNonFBAReactionsUB potrà dinamicamente aggiornare anche queste forward FBA.
		 */
		void loadAndApplyFBAReactionUpperBounds(
				std::vector<FBGLPK::LPprob>& vec_fluxb,
				const std::string& /*unused*/)
		{
				// 1) trova il file CSV
				const std::string path = firstExisting({
				    "ub_bounds_projected_gui.csv",
				    "ub_bounds_projected.csv",
				    "EX_upper_bounds_FBA.csv"
				});
				if (path.empty()) {
				    std::cerr << "[Warning][FBA] projected bounds CSV not found – UBs unchanged.\n";
				    return;
				}
				std::cout << "[DEBUG][FBA] Loading projected bounds from: " << path << "\n";

				std::ifstream file(path);
				if (!file.good()) {
				    std::cerr << "[Warning][FBA] Cannot open '" << path << "' – UBs unchanged.\n";
				    return;
				}

				using Key = std::pair<std::string, std::size_t>;
				using DebugRow = std::tuple<std::string, std::string, double, double, double>;
				std::unordered_set<Key, RxnLPHash> done;
				std::vector<DebugRow> debugRows;

				// 2) salta header
				std::string header;
				std::getline(file, header);
				std::cout << "[DEBUG][FBA] Skipped header: " << header << "\n";

				bool firstRow = !externalVolumeLoaded;
				std::string line;
				while (std::getline(file, line)) {
				    if (line.empty()) continue;
				    std::istringstream ss(line);

				    std::string reaction, model, concStr, volStr;
				    std::getline(ss, reaction, ',');
				    std::getline(ss, model,    ',');
				    std::getline(ss, concStr,  ',');
				    std::getline(ss, volStr,   ',');

				    if (reaction.empty() || concStr.empty() || volStr.empty()) {
				        std::cerr << "[Warning][FBA] Skipping invalid line: " << line << "\n";
				        continue;
				    }
				    std::cout << "[DEBUG][FBA] Parsed → reaction='" << reaction
				              << "', model='" << model
				              << "', conc='" << concStr
				              << "', vol='" << volStr << "'\n";

				    // 3) parse dei valori
				    double conc = 0.0, volume = 1.0;
				    try {
				        conc   = std::stod(concStr);
				        volume = std::stod(volStr);
				    } catch (...) {
				        std::cerr << "[Warning][FBA] non-numeric conc/volume for '" << reaction << "'\n";
				        continue;
				    }

				    // 4) cattura systemExternalVolume solo alla prima riga valida
				    if (firstRow) {
				        systemExternalVolume = volume;
				        externalVolumeLoaded = true;
				        firstRow = false;
				        std::cout << "[DEBUG][FBA] systemExternalVolume set to " << systemExternalVolume << "\n";
				    }

				    // 5) trova l’indice LP corrispondente
				    std::size_t lpIdx = findLPIndex(vec_fluxb, model);
				    if (lpIdx == std::numeric_limits<std::size_t>::max()) {
				        std::cerr << "[Warning][FBA] model '" << model
				                  << "' not found for reaction '" << reaction << "'\n";
				        continue;
				    }

				    Key key{reaction, lpIdx};
				    if (done.count(key)) {
				        std::cout << "[DEBUG][FBA] Already handled (" << reaction 
				                  << "," << vec_fluxb[lpIdx].getFilename() << ")\n";
				        continue;
				    }
				    done.insert(key);

				    // 6) consideriamo solo le forward (_f)
				    if (!endsWith(reaction, "_f")) {
				        std::cout << "[DEBUG][FBA] Skipping non-forward reaction: " << reaction << "\n";
				        continue;
				    }

				    // 7) detect “external” se prima di “_f” compare “_e”
				    bool isExternal = false;
				    auto pos = reaction.rfind("_f");
				    if (pos != std::string::npos && pos >= 2 &&
				        reaction.substr(pos-2,2) == "_e") {
				        isExternal = true;
				    }
				    std::string subtype = isExternal ? "exchange" : "internal";
				    std::cout << "[DEBUG][FBA] Reaction " << reaction
				              << (isExternal ? " is external\n" : " is internal\n");

				    // 8) applica subito il bound statico
				    int col = vec_fluxb[lpIdx].fromNametoid(reaction);
				    if (col < 0) {
				        std::cerr << "[Warning][FBA] reaction '" << reaction 
				                  << "' not in LP '" << vec_fluxb[lpIdx].getFilename() << "'\n";
				        continue;
				    }
				    double newUb = conc * systemExternalVolume;
				    double lb    = vec_fluxb[lpIdx].getLwBounds(col);
				    std::cout << "[DEBUG][FBA] Applying static bound on LP '"
				              << vec_fluxb[lpIdx].getFilename() << "' col=" << col
				              << ": lb=" << trunc(lb,decimalTrunc)
				              << ", ub=" << trunc(newUb,decimalTrunc) << "\n";
				    vec_fluxb[lpIdx].update_bound(
				        col,
				        "GLP_DB",
				        trunc(lb, decimalTrunc),
				        trunc(newUb, decimalTrunc)
				    );

				    // 9) registra in NonFBAReactionBaseUB & ReactionSubtype per update dinamico
				    NonFBAReactionBaseUB[key] = conc;
				    ReactionSubtype[key]      = subtype;
				    std::cout << "[DEBUG][FBA] Stored for dynamic update: key=("
				              << reaction << ",LP#" << lpIdx
				              << "), C_m=" << conc
				              << ", subtype=" << subtype << "\n";

				    // 10) colleziona per CSV
				    debugRows.emplace_back(
				        reaction,
				        vec_fluxb[lpIdx].getFilename(),
				        conc,
				        volume,
				        newUb
				    );
				}
				file.close();

				// 11) registra projectedFBA per eventuali filtri
				projectedFBA.clear();
				for (auto& tup : debugRows) {
				    projectedFBA.insert(std::get<0>(tup));
				}
				std::cout << "[DEBUG][FBA] Collected " << projectedFBA.size()
				          << " projectedFBA reactions\n";

				// 12) dump CSV di debug
				/*if (!debugRows.empty()) {
				    std::ofstream dbg("debug_FBA_projected_bounds.csv");
				    dbg << "reaction,LPfile,conc,volume,new_upper_bound\n";
				    for (auto &r : debugRows) {
				        dbg << std::get<0>(r) << ','
				            << std::get<1>(r) << ','
				            << std::get<2>(r) << ','
				            << std::get<3>(r) << ','
				            << std::get<4>(r) << '\n';
				    }
				    dbg.close();
				    std::cout << "[DEBUG][FBA] Wrote debug_FBA_projected_bounds.csv\n";
				} else {
				    std::cout << "[DEBUG][FBA] No projected bounds to write\n";
				}*/
		}


		/*─────────────────────────────────────────────────────────────────────────────
		 *  updateNonFBAReactionsUB
		 *
		 *  - gestisce **solo** le reazioni NON-projected (_r e _f)
		 *  - per le _f mantiene l’UB fisso letto dal CSV
		 *  - per le _r calcola l’UB dinamico:
		 *
		 *         newUB = |Cₘ| · λ
		 *
		 *         λ =  1 / Σ_k (xN_k · xB_k · 10⁻¹²)      se subtype == "exchange"
		 *             1 /   (xN_s · xB_s · 10⁻¹²)         altrimenti (demand/sink int.)
		 *
		 *  dove Cₘ è il valore (concentrazione) salvato in NonFBAReactionBaseUB.
		 *  MW_B = 1 quindi il fattore 10⁻¹² converte pgDW → gDW.
		 *────────────────────────────────────────────────────────────────────────────*/
			void updateNonFBAReactionsUB(
					std::vector<FBGLPK::LPprob>& vec_fluxb,
					std::map<std::string,int>&   NumPlaces,
					double*                      Value)
			{
					using Row = std::tuple<std::string,std::string,double,double,double,double,double>;
					std::vector<Row> dbg;

				//	std::cout << "[DEBUG updateNonFBA] Starting dynamic update of non-projected bounds (_f & _r)\n";

					auto lpListFor = [&](const std::string& rxn){
						  std::vector<std::size_t> v;
						  for (auto& kv : NonFBAReactionBaseUB)
						      if (kv.first.first == rxn)
						          v.push_back(kv.first.second);
						  if (!v.empty()) return v;
						  if (reactionToFileMap.count(rxn))
						      for (auto i : reactionToFileMap[rxn])
						          v.push_back(i);
						  return v;
					};

					std::unordered_set<std::string> allRxn;
					for (auto& kv : reactionToFileMap)    allRxn.insert(kv.first);
					for (auto& kv : NonFBAReactionBaseUB) allRxn.insert(kv.first.first);

					for (const auto& rxn : allRxn) {
						  // 1) skip solo le reverse FBA (_r) — forward FBA (_f) **vengono** elaborate qui
						  if (FBAreactions.count(rxn) && endsWith(rxn, "_r")) {
						     // std::cout << "[DEBUG updateNonFBA] Skipping reverse FBA-bound: " << rxn << "\n";
						      continue;
						  }
						  // 2) skip le biomassa non-FBA / FBA
						  if (rxn == "EX_biomass_e_f" || rxn == "EX_biomass_e_r") {
						    //  std::cout << "[DEBUG updateNonFBA] Skipping biomass reaction: " << rxn << "\n";
						      continue;
						  }

						//  std::cout << "[DEBUG updateNonFBA] Processing reaction: " << rxn << "\n";
						  auto lpIdxList = lpListFor(rxn);
						  if (lpIdxList.empty()) {
						  //    std::cout << "[DEBUG updateNonFBA]   No LPs found for " << rxn << "\n";
						      continue;
						  }

						  // prendi baseUb (C_m) e tipo (exchange/internal) dal primo LP
						  auto key0   = std::make_pair(rxn, lpIdxList[0]);
						  double C_m  = fabs(NonFBAReactionBaseUB[key0]);
						  std::string subtype = ReactionSubtype.count(key0)
						                        ? ReactionSubtype[key0]
						                        : "exchange";
						  //std::cout << "[DEBUG updateNonFBA]   baseUb C_m=" << C_m
						   //         << ", subtype=" << subtype << "\n";

						  // calcola λ sommando su tutti i LP (exchange) o solo sul primo (internal)
						  double denom = 0.0;
						  if (subtype == "exchange") {
						     // std::cout << "[DEBUG updateNonFBA]   exchange: sum pop·biomass\n";
						      for (auto idx : lpIdxList) {
						          double pop = problemBacteriaPlace.count(idx)
						                       ? std::floor(Value[ NumPlaces.at(problemBacteriaPlace[idx]) ])
						                       : 1.0;
						          double biomass = problemBiomassPlace.count(idx)
						                       ? trunc(Value[ NumPlaces.at(problemBiomassPlace[idx]) ], decimalTrunc)
						                       : 1.0;
						          double contrib = pop * biomass * 1e-12;
						       /*   std::cout << "[DEBUG updateNonFBA]     LP#" << idx
						                    << ": pop=" << pop
						                    << ", biomass=" << biomass
						                    << " → contrib=" << contrib << "\n";*/
						          denom += contrib;
						      }
						  } else {
						      size_t idx = lpIdxList[0];
						      double pop = problemBacteriaPlace.count(idx)
						                   ? std::floor(Value[ NumPlaces.at(problemBacteriaPlace[idx]) ])
						                   : 1.0;
						      double biomass = problemBiomassPlace.count(idx)
						                   ? trunc(Value[ NumPlaces.at(problemBiomassPlace[idx]) ], decimalTrunc)
						                   : 1.0;
						      denom = pop * biomass * 1e-12;
						     /* std::cout << "[DEBUG updateNonFBA]   internal LP#" << idx
						                << ": pop=" << pop
						                << ", biomass=" << biomass
						                << " → denom=" << denom << "\n";*/
						  }

						  double lambda = denom > 0.0 ? 1.0/denom : 0.0;
						  //std::cout << "[DEBUG updateNonFBA]   Total denom=" << denom
						  //          << " → lambda=" << lambda << "\n";

						  // applichiamo il bound dinamico sia alle forward (_f) che alle reverse (_r) non-FBA
						  // e alle forward FBA (_f), perché qui abbiamo tolto il filtro su projectFBA _f
						  for (auto idx : lpIdxList) {
						      int col = vec_fluxb[idx].fromNametoid(rxn);
						      if (col < 0) continue;
						      double oldUb = vec_fluxb[idx].getUpBounds(col);
						      double newUb = C_m * lambda;
						     /* std::cout << "[DEBUG updateNonFBA]   LP#" << idx
						                << ": oldUb=" << oldUb
						                << ", newUb=" << newUb << "\n";*/
						      if (newUb <= 0.0) {
						         // std::cout << "[DEBUG updateNonFBA]     Setting GLP_FX 0\n";
						          vec_fluxb[idx].update_bound(col, "GLP_FX", 0.0, 0.0);
						      } else {
						         // std::cout << "[DEBUG updateNonFBA]     Setting GLP_DB [0," 
						          //          << trunc(newUb,decimalTrunc) << "]\n";
						          vec_fluxb[idx].update_bound(col, "GLP_DB", 0.0, trunc(newUb,decimalTrunc));
						      }
						      dbg.emplace_back(
						          rxn,
						          vec_fluxb[idx].getFilename(),
						          problemBacteriaPlace.count(idx)
						              ? std::floor(Value[ NumPlaces.at(problemBacteriaPlace[idx]) ])
						              : 1.0,
						          problemBiomassPlace.count(idx)
						              ? trunc(Value[ NumPlaces.at(problemBiomassPlace[idx]) ], decimalTrunc)
						              : 1.0,
						          C_m,
						          oldUb,
						          newUb
						      );
						  }
					}

				/*	std::cout << "[DEBUG updateNonFBA] Writing CSV 'debug_nonFBA_runtime_bounds.csv'\n";
					std::ofstream out("debug_nonFBA_runtime_bounds.csv");
					out << "reaction,LPfile,pop,biomass,Cm(oldBase),oldUB,newUB\n";
					for (auto& r : dbg) {
						  out << std::get<0>(r) << ','
						      << std::get<1>(r) << ','
						      << std::get<2>(r) << ','
						      << std::get<3>(r) << ','
						      << std::get<4>(r) << ','
						      << std::get<5>(r) << ','
						      << std::get<6>(r) << '\n';
					}*/
					//out.close();
					//std::cout << "[DEBUG updateNonFBA] Dynamic update complete\n";
			}


		void updateAllBiomassReactionsUpperBounds(
				    double*                        Value,
				    const map<string,int>&         NumPlaces,
				    vector<FBGLPK::LPprob>&        vec_fluxb)
		{
				problemsWithLowBiomass.clear();

				for (const auto& transition : biomassTransitions) {

				    const size_t problemIndex = FBAproblems[transition];
				    const double BioMax       = vec_fluxb[problemIndex].getBioMax();   // x_B,max [pgDW/cell]
				    const double BioMin       = vec_fluxb[problemIndex].getBioMin();   // x_B,min [pgDW/cell]
				    const double xB           = trunc(Value[ NumPlaces.at( problemBiomassPlace.at(problemIndex) ) ],
				                                      decimalTrunc);                   // x_B    [pgDW/cell]

				    const double Gamma        = BioMax - xB;                           // Γ      [pgDW/cell]
				    const double psi        = 1.0;                           					 // psi     [pgDW/cell]

				    /* ---------- calcolo nuovo UB secondo Eq.(5.16) ---------- */
					double newUB = 0.0;
					double mu = muMaxMap.count(problemIndex) ? muMaxMap[problemIndex] : 1.0;

					if (xB < BioMin) {
							newUB = 0.0;                                  // metabolic distress
							problemsWithLowBiomass[problemIndex] = true;
					}
					else if (Gamma > 0.0) {
							newUB = (Gamma / psi)  * mu;                  			// mmol / gDW / h
					}
					else {                                            // Γ ≤ 0
							newUB = Lcutoff;                              // growth-limitation cut-off
					}

				    /* ---------- applicazione bound ---------- */
				    const int rxnIdx   = vec_fluxb[problemIndex].getPFBA_index();  // indice reazione biomassa
				    const double LBold = vec_fluxb[problemIndex].getLwBounds(rxnIdx);

				    if (newUB <= LBold + 1e-15) {                  // ≈ stesso valore → fixa a 0
				        vec_fluxb[problemIndex].update_bound(rxnIdx, "GLP_FX", 0.0, 0.0);
				    } else {
				        vec_fluxb[problemIndex].update_bound(rxnIdx, "GLP_DB",
				                                             LBold,
				                                             trunc(newUB, decimalTrunc));
				    }

				    /* ---------- DEBUG verbose ---------- */
						std::cout << "[Biomass UB] LP="  << vec_fluxb[problemIndex].getFilename()
											<< "  xB="             << xB
											<< "  BioMin="         << BioMin
											<< "  BioMax="         << BioMax
											<< "  Γ="              << Gamma
											<< "  μ_max="          << mu
											<< "  → newUB="        << newUB << " mmol/gDW/h\n";

				}
		}


    void performPFBA(vector<class FBGLPK::LPprob>& vec_fluxb, int problemIndex) {
        double optimalBiomass = vec_fluxb[problemIndex].getOBJ();
        
        int biomassIndex = vec_fluxb[problemIndex].getPFBA_index();
        double originalLb = vec_fluxb[problemIndex].getLwBounds(biomassIndex);
        double originalUb = vec_fluxb[problemIndex].getUpBounds(biomassIndex);
        int originalType = vec_fluxb[problemIndex].get_bound_type(biomassIndex);
        
        // Prove debug metto un piccolo range
        double eps = 1e-6;
        double lb = optimalBiomass - eps;
        double ub = optimalBiomass + eps;
        vec_fluxb[problemIndex].update_bound(biomassIndex, "GLP_DB", lb, ub);
        
        vec_fluxb[problemIndex].setMinimizeFluxObjective(biomassIndex);
        
        // Rieseguo il solver con la pFBA
        //cout << "Ora risolvo pFBA" << endl;
        vec_fluxb[problemIndex].solve();
        //cout << "RISOLTO pFBA" << endl;
        
        // -- DEBUG 3: Stampo la situazione pFBA
        //vec_fluxb[problemIndex].debugPFBA();
        
        // Step 3: Ripristino i bound originali
        vec_fluxb[problemIndex].update_bound(biomassIndex, originalType, originalLb, originalUb);
        
        // Step 3.1: Torno a massimizzare la biomassa
        vec_fluxb[problemIndex].resetMaximizationObjectiveForBiomass(biomassIndex);
    }

		/*───────────────────────────────────────────────────────────────
		 *  loadMuMaxValues
		 *  - CSV con header:  Model,mu_max
		 *  - se non esiste o vuoto → μmax = 1 per tutti i modelli
		 *───────────────────────────────────────────────────────────────*/
		void loadMuMaxValues(const std::vector<FBGLPK::LPprob>& vec_fluxb)
		{
				const std::string path = firstExisting({
				    "mu_max_values_gui.csv",
				    "mu_max_values.csv"
				});

				/* default: 1 per tutti i modelli */
				muMaxMap.clear();
				for (std::size_t i = 0; i < vec_fluxb.size(); ++i)
				    muMaxMap[i] = 1.0;

				if (path.empty()) {
				    std::cerr << "[μmax] CSV not found – using μmax = 1 for every model.\n";
				    return;
				}

				std::ifstream file(path.c_str());
				if (!file.good()) {
				    std::cerr << "[μmax] Cannot open file '" << path
				              << "' – using μmax = 1 for every model.\n";
				    return;
				}

				std::string header;
				std::getline(file, header);             // salta header

				std::string line;
				std::size_t valid = 0;
				while (std::getline(file, line)) {
				    if (line.empty()) continue;

				    std::istringstream ss(line);
				    std::string model, muStr;
				    std::getline(ss, model, ',');
				    std::getline(ss, muStr, ',');

				    if (model.empty() || muStr.empty()) continue;

				    double mu;
				    try { mu = std::stod(muStr); }
				    catch (...) {
				        std::cerr << "[μmax] Warning: mu_max non numerico per modello '"
				                  << model << "'\n";
				        continue;
				    }

				    std::size_t idx = findLPIndex(vec_fluxb, model);
				    if (idx == std::numeric_limits<std::size_t>::max()) {
				        std::cerr << "[μmax] Warning: modello '" << model
				                  << "' non trovato fra i LP – skip\n";
				        continue;
				    }

				    muMaxMap[idx] = mu;
				    ++valid;
				}
				file.close();

				if (valid == 0) {
				    std::cerr << "[μmax] File '" << path
				              << "' vuoto o senza righe valide – μmax = 1 per tutti i modelli.\n";
				} else {
				    std::cout << "[μmax] Loaded " << valid
				              << " μmax value(s) from '" << path << "'.\n";
				}
		}
		
		// 1) Debug‐enhanced mapMetabolitesToReactions()
		//    Call this right after you’ve built FBAplace (in init_data_structures_class).
		void mapMetabolitesToReactions() {
				metaboliteToReactions.clear();
				std::cout << "[DEBUG mapMetabolitesToReactions] FBAplace size = " 
				          << FBAplace.size() << "\n";

				for (const auto &kv : FBAplace) {
				    const std::string &reaction  = kv.first;
				    const std::string &placesStr = kv.second;
				    auto places = splitAndTrim(placesStr, ',');

				    std::cout << "[DEBUG mapMetabolitesToReactions] reaction = " 
				              << reaction << " → places: ";
				    for (auto &p : places) std::cout << p << ", ";
				    std::cout << "\n";

				    for (auto &p : places) {
				        metaboliteToReactions[p].insert(reaction);
				    }
				}

				// print summary
				std::cout << "[DEBUG metaboliteToReactions] summary:\n";
				for (const auto &kv : metaboliteToReactions) {
				    std::cout << "  metabolite " << kv.first << " → reactions: ";
				    for (auto &r : kv.second) std::cout << r << ", ";
				    std::cout << "\n";
				}
		}



};
