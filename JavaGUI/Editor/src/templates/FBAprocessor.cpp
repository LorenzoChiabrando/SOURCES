
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
    double process(double *Value, vector<class FBGLPK::LPprob>& vec_fluxb, map<string,int>& NumTrans, map<string,int>& NumPlaces, const vector<string>& NameTrans, const struct InfTr* Trans, int T, const double& time) {
        if (!init) {
            init_data_structures_class(NameTrans, vec_fluxb, Value, NumPlaces); // Ensure it is initialized before proceeding
            mapReactionsFromProblems(vec_fluxb);
            loadAndApplyFBAReactionUpperBounds(vec_fluxb, "EX_upper_bounds_FBA.csv");
            updateNonFBAReactionUpperBoundsFromFile(vec_fluxb, "EX_upper_bounds_nonFBA.csv");
            loadMuMaxValues(vec_fluxb);
            loadGeneRules(vec_fluxb, "GeneRules.txt");
            debugPrintGeneRules();
            firstTransitionName = NameTrans[T];
            init = true;
        }
        if (firstTransitionName == NameTrans[T]) {
            for (auto& problem : FBAproblems) {
                if (hasMultiSpecies && floor(Value[NumPlaces.at(problemBacteriaPlace.at(problem.second))]) < 1) {
                    deadBacterialSpecies[problem.second] = true;
                }
            }
            if(hasBioMASS){
            		updateAllBiomassReactionsUpperBounds(Value, NumPlaces, vec_fluxb); // Update biomass upper limits
            }
            updateFluxBoundsAndSolve(Value, vec_fluxb, NumPlaces, time); // Update fluxes only on first transition
        }
        double rate = 0;
        size_t problemIndex = FBAproblems[NameTrans[T]];
        // Return zero rate if the species for this transition is dead
        if (deadBacterialSpecies.find(problemIndex) != deadBacterialSpecies.end()) {
            return 0; // Return zero rate for dead species
        }

        // Check for transitions without associated file or specific biomass transitions.
        if (transitionsWithoutFile.find(NameTrans[T]) != transitionsWithoutFile.end()) {
            return rate; // Skip further processing for these cases
        }

        // Default computation if none of the above conditions met.
        rate = computeRate(vec_fluxb, NumPlaces, NameTrans, Value, T, decimalTrunc, time); // Compute and return the rate for the given transition
        return rate; // Return the calculated rate
    }



		private:
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

						string finalReaction;
            if (!lpFile.empty() && !reaction.empty()) {
								 size_t index = findLPIndex(vec_fluxb, lpFile);
								 if (index == numeric_limits<size_t>::max()) {
    						 		// Non trovata => transitionsWithoutFile
    								transitionsWithoutFile.insert(transition);
								 } else {
    								const auto& forwardMap  = vec_fluxb[index].getForwardReactions();
    								const auto& reverseMap  = vec_fluxb[index].getReverseReactions();
    								const auto& irreversMap = vec_fluxb[index].getIrreversibileReactions();

										//cerr << "[DEBUG TRANSIZIONE] '" << transition << endl;
										  finalReaction = standardizeReactionName(
    									transition,      
    									reaction,        
    									inputs,
    									outputs,
    									forwardMap,
    									reverseMap,
    									irreversMap,
   	 								isBiomass
										);
    								//cout << "[DEBUG readFBAInfo] reactionBase='" << reaction
         				  //<< "transition: " << transition	<< "' => finalReaction='" << finalReaction << "'\n";

    								FBAreact[transition] = finalReaction;
                   int checkReacId = vec_fluxb[index].fromNametoid(finalReaction);
                   if (checkReacId == -1) {
                   		cerr << "[WARNING readFBAInfo] Reaction '" << finalReaction
                            << "' was NOT found in LP index=" << index << " (file='" << lpFile << "').\n";
        								return false;
                   }

    							if (!bacteriaCountPlace.empty() && bacteriaCountPlace != "N/A") {
        							problemBacteriaPlace[index] = bacteriaCountPlace;
    							}
    							if (!bacteriaBiomassPlace.empty() && bacteriaBiomassPlace != "N/A") {
        							problemBiomassPlace[index] = bacteriaBiomassPlace;
    							}

    							problems.insert(index);
    							bacteriaToBioMax[index]  = vec_fluxb[index].getBioMax();
    							bacteriaToBioMean[index] = vec_fluxb[index].getBioMean();
    							bacteriaToBioMin[index]  = vec_fluxb[index].getBioMin();

    							FBAproblems[transition] = index;
    							problemsToReactions[index].insert(transition);
								 }
            		 if (isBiomass) {
                		FBAreactions.insert(reaction + "_f");
                		FBAreactions.insert(reaction + "_r");
                }else{
                		FBAreactions.insert(finalReaction);
                }
                if (FBAreact.find(transition) != FBAreact.end()) {
                    string reactionName = FBAreact[transition];
                    if ((!inputs.empty() || !outputs.empty())) {
                        string combinedPlaces;
                        for (const auto& place : inputs) {
                            if (!place.empty() && combinedPlaces.find(place + ",") == string::npos) {
                                combinedPlaces += place + ",";
                            }
                        }
                        for (const auto& place : outputs) {
                            if (!place.empty() && combinedPlaces.find(place + ",") == string::npos) {
                                combinedPlaces += place + ",";
                            }
                        }
                        if (!combinedPlaces.empty()) {
                            combinedPlaces.pop_back(); // Remove the last comma
                            FBAplace[reactionName] = combinedPlaces;
                        }
                    }
                }
            }
            ReactionMultiplicity[transition] = multiplicity;
            //cout << "For transition: " << transition << " multiplicity: " << ReactionMultiplicity[transition] << endl; 

            // Populate the reactionToBacteria map if bacteriaCountPlace is specified
            if (!bacteriaCountPlace.empty() && bacteriaCountPlace != "N/A") {
                reactionToBacteria[transition] = bacteriaCountPlace;
                hasMultiSpecies = true;
            }
            if (!bacteriaBiomassPlace.empty() && bacteriaBiomassPlace != "N/A") {
                reactionToBacteriaBIOMASS[transition] = bacteriaBiomassPlace;
                hasBioMASS = true;
            }
        }

        return true;
    }


		double dynamicThreshold(double previousConcentration) {
    		double baseThreshold = 1e-3;
    		return std::max(minDeltaThreshold, baseThreshold * fabs(previousConcentration));
    		// return minDeltaThreshold;
		}


    set<string> checkSignificantChange(double* Value, map<string, int>& NumPlaces, double time) {
        reactionsToUpdate.clear();
        set<string> changedMetabolites;
        double relChange = (relEpsilon > 100.0) ? 100.0 : ((relEpsilon < 0.0 && absEpsilon == -1) ? 0 : ((relEpsilon < 0.0 && absEpsilon != -1) ? -1 : relEpsilon));
       //cout << "[DEBUG] relChange: " << relChange << endl;
       // cout << "[DEBUG] absChange: " << absEpsilon << endl;
        bool useAbsolute = (absEpsilon != -1);
        bool useRelative = (relChange != -1);
        for (const auto& place : FBAplace) {
            string metabolite = place.second;
            double currentConcentration = trunc(Value[NumPlaces[metabolite]], decimalTrunc);
            double previousConcentration = previousConcentrations[metabolite];
            double absoluteDiff = fabs(currentConcentration - previousConcentration);
            double percentChange = (previousConcentration != 0)
                ? (absoluteDiff / fabs(previousConcentration)) * 100.0
                : 100.0;
            bool condAbsolute = useAbsolute && (absoluteDiff > absEpsilon);
            bool condRelative = useRelative && (percentChange > relChange);
            bool condFromZero = useRelative && (previousConcentration == 0 && currentConcentration != 0);
            bool condFirstSol = (previousConcentration == -1);
            if (condAbsolute || condRelative || condFromZero || condFirstSol) {
                /*cout << "[DEBUG] Metabolite '" << metabolite
                     << "' changed. Old: " << previousConcentration
                     << ", New: " << currentConcentration
                     << ", Diff: " << (currentConcentration - previousConcentration)
                     << ", AbsDiff: " << absoluteDiff
                     << ", PercentChange: " << percentChange << "%"
                     << ", condAbsolute: " << condAbsolute
                     << ", condRelative: " << condRelative
                     << ", condFromZero: " << condFromZero << endl;*/
                changedMetabolites.insert(metabolite);
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
        double *Value, 
        vector<class FBGLPK::LPprob>& vec_fluxb, 
        map<string, int>& NumPlaces, 
        set<string>& changedMetabolites,
        double time
    ) {
        set<size_t> problemsToUpdate;

        // (1) Costruiamo l'insieme dei problemIndex da aggiornare
        for (const string& metabolite : changedMetabolites) {
            if (metaboliteToProblems.find(metabolite) != metaboliteToProblems.end()) {
                for (auto prob : metaboliteToProblems[metabolite]) {
                    if (deadBacterialSpecies.find(prob) == deadBacterialSpecies.end() 
                        || !deadBacterialSpecies[prob]) {
                        problemsToUpdate.insert(prob);
                    }
                }
            }
        }

         for (auto problemIndex : problemsToUpdate)
						for (const std::string& reaction : reactionsToUpdate)
						{
								/* ---- filtri invariati ---- */
								int colIdx = vec_fluxb[problemIndex].fromNametoid(reaction);
								if (colIdx == -1)                   continue;
								if (reaction == "EX_biomass_e_f" ||
								    reaction == "EX_biomass_e_r")  continue;
								if (!endsWith(reaction,"_r"))      continue;   // solo uptake

								/* ---- ①  pool disponibile  ------------------------------- */
								double conc_mmol_per_mL =
								    trunc( Value[ NumPlaces.at( FBAplace[reaction] ) ],
								           decimalTrunc);                      // C_m

								/* ---- ②  λ (collettivo o singolo) ------------------------ */
								calculateMultiplicativeConstant(Value, NumPlaces,
								                                problemIndex, reaction, vec_fluxb);
								double lambda = multiplicativeConstant;        // mL / (gDW·h)

								/* ---- ③  flusso orario teorico --------------------------- */
								double vr_max_h = conc_mmol_per_mL * lambda;   // mmol / (gDW·h)

								/* ---- ④  scala al passo d’integrazione ------------------ */
								double vr_max_dt = vr_max_h * INTEGRATION_STEP_H;   // mmol / (gDW)

								/* ---- ⑤  protegge dal pool-negativo --------------------- */
								if (conc_mmol_per_mL <= 0.0 || lambda == 0.0)
								    vr_max_dt = 0.0;

								/* ---- ⑥  applica il nuovo UB ---------------------------- */
								vec_fluxb[problemIndex].update_bound(
								    colIdx,
								    (vr_max_dt <= 0.0 ? "GLP_FX" : "GLP_DB"),
								    0.0,
								    trunc(vr_max_dt, decimalTrunc)
								);

								std::cout << "[UB-upd] " << reaction
								          << "  LP="   << vec_fluxb[problemIndex].getFilename()
								          << "  C_m="  << conc_mmol_per_mL
								          << "  λ="    << lambda
								          << "  UB(h)="<< vr_max_h
								          << "  UB(dt)="<< vr_max_dt << '\n';
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
        }
    }



    /**
     * Updates the flux bounds and solves the FBA problems only for those metabolites that have undergone significant changes.
     * This method first checks for significant changes in metabolite concentrations. If changes are detected,
     * it updates the flux bounds for those metabolites, solves the corresponding FBA problems,
     * and then updates the metabolite concentrations in the model.
     *
     * @param Value Pointer to an array containing current metabolic values.
     * @param vec_fluxb Reference to a vector of FBGLPK::LPprob objects, each representing a distinct FBA problem.
     * @param NumPlaces Mapping of place names to their indices in the metabolic array.
     * @param time Current time to manage updates efficiently.
     */
    void updateFluxBoundsAndSolve(double* Value, vector<class FBGLPK::LPprob>& vec_fluxb, map<string, int>& NumPlaces, double time) {
        set<string> changedMetabolites = checkSignificantChange(Value, NumPlaces, time);
        if (!changedMetabolites.empty()) {
            for (const auto& place : FBAplace) {
                reactionsToUpdate.insert(place.first);
            }

            applyGeneRegulationRules(Value, vec_fluxb, NumPlaces, time);
            updateFluxBounds(Value, vec_fluxb, NumPlaces, changedMetabolites, time);
            updateNonFBAReactionsUB(vec_fluxb, NumPlaces, Value);
            solveFBAProblems(vec_fluxb, changedMetabolites);
            updateConcentrations(Value, NumPlaces, changedMetabolites);
            count += 1;
            //cout << "risolvo per la: " << count << endl;
        }
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

					/*―― 1. ricavo il flusso specifico ottimo v★ [mmol/gDW/h] ――*/
					int varIndex   = vec_fluxb[problemIndex].fromNametoid( FBAreact[transitionName] );
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
					if (isBiomass) {
						  /*  ζ_B  =  MW_B · x_B   (MW_B = 1 gDW/mmol)                       */
						  zeta = xB;                                    // pgDW / cell
					} else {
						  /*  ζ(m) = x_N · x_B · 10⁻¹²  [Eq.(37)]                            */
						  zeta = static_cast<double>(xN) * xB * 1e-12;  // gDW (⇒ v★·ζ → mmol/h)
					}

					/*―― 4. rate finale ――*/
					double rate = vStar * zeta * mult;

					/*―― 5. DEBUG VERBOSO ――*/
					std::cout << "\n[computeRate] transition=" << transitionName
						        << (isBiomass ? "  (biomass)" : "  (boundary)")
						        << "\n   v*                = " << vStar               << "  [mmol/gDW/h]"
						        << "\n   x_B               = " << xB                  << "  [pgDW/cell]"
						        << "\n   x_N               = " << xN                  << "  [cell]"
						        << "\n   ζ                 = " << zeta
						        << (isBiomass ? "  [pgDW/cell]" : "  [gDW]")
						        << "\n   multiplicity (μ)  = " << mult
						        << "\n   → rate            = " << rate
						        << (isBiomass ? "  [pgDW/cell/h]" : "  [mmol/h]")
						        << std::endl;

					return rate;
			}



    /**
     * Initializes the metabolic concentrations from provided values.
     * This method initializes the concentrations for each metabolite based on the indices provided in a map.
     * Each concentration is truncated to a specified number of decimal places before being set.
     *
     * @param Value Array of initial concentration values.
     * @param NumPlaces Mapping from metabolite names to their indices in the Value array.
     * @return void
     */
    void initializeConcentrations(double* Value, map<string, int>& NumPlaces) {
        for (const auto& pair : FBAplace) {
            string placeName = pair.second;
            previousConcentrations[pair.second] = -1;
        }
    }


    /**
     * Maps each metabolite to the set of LP problems that it affects.
     * This method initializes the mapping from metabolites to their related LP problems,
     * which is used to optimize updates when only certain metabolites change.
     */
    void mapMetabolitesToProblems() {
        for (const auto& problem : FBAproblems) {
            const string& reaction = FBAreact[problem.first];
            if (FBAplace.find(reaction) != FBAplace.end()) {
                size_t problemIndex = problem.second;
                string places = FBAplace[reaction];
                vector<string> metabolites = splitAndTrim(places, ',');
                for (const string& metabolite : metabolites) {
                    metaboliteToProblems[metabolite].insert(problemIndex);
                }
            }
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

		/*─────────────────────────────────────────────────────────────────────────────
		 *  Lettura CSV upper-bounds NON-projected
		 *  – memorizza UB di base  (NonFBAReactionBaseUB)
		 *  – salva subtype         (ReactionSubtype)   ←  "exchange" | "internal"
		 *────────────────────────────────────────────────────────────────────────────*/
		void updateNonFBAReactionUpperBoundsFromFile(
				    std::vector<FBGLPK::LPprob>& vec_fluxb,
				    const std::string& /*unused*/)
		{
				/*─────────────── percorsi prioritari ───────────────────────────*/
				const std::string pathFwd = firstExisting({
				    "non_projected_forward_bounds_gui.csv",
				    "non_projected_forward_bounds.csv",
				    "ub_bounds_not_projected_gui.csv",      // legacy fallback
				    "ub_bounds_not_projected.csv",
				    "EX_upper_bounds_nonFBA.csv"            // last‑chance fallback
				});
				const std::string pathRev = firstExisting({
				    "non_projected_reverse_background_met_gui.csv",
				    "non_projected_reverse_background_met.csv"
				});

				if (pathFwd.empty() && pathRev.empty()) {
				    std::cerr << "Warning: non‑FBA bounds CSVs not found – UBs unchanged." << std::endl;
				    return;                                      // niente da fare
				}

				using Key  = std::pair<std::string,std::size_t>;         // (reaction,LPidx)
				using RowF = std::tuple<std::string,std::string,double>; // debug riga forward
				using RowR = std::tuple<std::string,std::string,double,double,double>; // reverse

				std::unordered_set<Key,RxnLPHash> done;                  // evita duplicati
				std::vector<RowF> debugFwd;
				std::vector<RowR> debugRev;

				auto parseLpList = [&](const std::string& rxn,
				                       const std::string& model)->std::vector<std::size_t>
				{
				    std::vector<std::size_t> lpList;
				    /* 1) match esplicito */
				    if (!model.empty()) {
				        size_t idx = findLPIndex(vec_fluxb, model);
				        if (idx != std::numeric_limits<size_t>::max()) lpList.push_back(idx);
				    }
				    /* 2) fallback reaction→LP map */
				    if (lpList.empty()) {
				        auto it = reactionToFileMap.find(rxn);
				        if (it != reactionToFileMap.end())
				            lpList.insert(lpList.end(), it->second.begin(), it->second.end());
				    }
				    return lpList;
				};

				/*──────────────────────── 1) CSV FORWARD (_f) ───────────────────────*/
				if (!pathFwd.empty()) {
				    std::ifstream file(pathFwd.c_str());
				    if (!file.good()) {
				        std::cerr << "[NonFBA‑FWD] Impossibile aprire '"<<pathFwd<<"' – skip"<<std::endl;
				    } else {
				        std::string header;  std::getline(file, header);   // header ignorato
				        std::string line;
				        while (std::getline(file,line)) {
				            if (line.empty()) continue;
				            std::istringstream ss(line);
				            std::string reaction, model, ubStr;
				            std::getline(ss,reaction,',');
				            std::getline(ss,model,',');
				            std::getline(ss,ubStr,',');
				            if (reaction.empty() || ubStr.empty()) continue;

				            double baseUb = 0.0;
				            try { baseUb = std::stod(ubStr); }
				            catch (...) { std::cerr << "[NonFBA‑FWD] UB non numerico per "<<reaction<<"\n"; continue; }

				            const auto lpList = parseLpList(reaction,model);
				            if (lpList.empty()) {
				                std::cerr << "[NonFBA‑FWD] nessun LP trovato per '"<<reaction<<"'\n";
				                continue;
				            }

				            const std::string subtype = (reaction.rfind("EX_",0)==0) ? "exchange" : "internal";
				            for (size_t idx: lpList) {
				                Key key{reaction,idx};
				                if (done.count(key)) continue;
				                done.insert(key);

				                NonFBAReactionBaseUB[key] = baseUb;
				                ReactionSubtype[key]      = subtype;
				                debugFwd.emplace_back(reaction, vec_fluxb[idx].getFilename(), baseUb);
				            }
				        }
				        file.close();
				    }
				}

				/*──────────────────────── 2) CSV REVERSE (_r) ───────────────────────*/
				if (!pathRev.empty()) {
				    std::ifstream file(pathRev.c_str());
				    if (!file.good()) {
				        std::cerr << "[NonFBA‑REV] Impossibile aprire '"<<pathRev<<"' – skip"<<std::endl;
				    } else {
				        std::string header;  std::getline(file, header);
				        std::string line;
				        while (std::getline(file,line)) {
				            if (line.empty()) continue;
				            std::istringstream ss(line);
				            std::string reaction, model, bgStr, volStr;
				            std::getline(ss,reaction,',');
				            std::getline(ss,model,',');
				            std::getline(ss,bgStr,',');
				            std::getline(ss,volStr,',');
				            if (reaction.empty() || bgStr.empty() || volStr.empty()) continue;

				            double background=0.0, volume=0.0;
				            try {
				                background = std::stod(bgStr);
				                volume     = std::stod(volStr);
				            } catch(...) {
				                std::cerr << "[NonFBA‑REV] valori non numerici per "<<reaction<<"\n";
				                continue;
				            }
				            double baseUb = background * volume;          // mmol

				            const auto lpList = parseLpList(reaction,model);
				            if (lpList.empty()) {
				                std::cerr << "[NonFBA‑REV] nessun LP trovato per '"<<reaction<<"'\n";
				                continue;
				            }
				            const std::string subtype = (reaction.rfind("EX_",0)==0) ? "exchange" : "internal";
				            for (size_t idx: lpList) {
				                Key key{reaction,idx};
				                if (done.count(key)) continue;
				                done.insert(key);

				                NonFBAReactionBaseUB[key] = baseUb;
				                ReactionSubtype[key]      = subtype;
				                debugRev.emplace_back(reaction, vec_fluxb[idx].getFilename(), background, volume, baseUb);
				            }
				        }
				        file.close();
				    }
				}

				/*──────────────────── nessuna riga valida? → early exit ─────────────*/
				if (debugFwd.empty() && debugRev.empty()) {
				    std::cerr << "[NonFBA] CSV presenti ma nessuna riga valida – UBs inalterati." << std::endl;
				    return;
				}

				/*──────────────────────── dump debug CSV ─────────────────────────────*/
				if (!debugFwd.empty()) {
				    std::ofstream dbg("debug_nonFBA_forward_bounds.csv");
				    dbg << "reaction,LPfile,baseUB\n";
				    for (auto& r: debugFwd)
				        dbg << std::get<0>(r) << ',' << std::get<1>(r) << ',' << std::get<2>(r) << '\n';
				    dbg.close();
				}
				if (!debugRev.empty()) {
				    std::ofstream dbg("debug_nonFBA_reverse_bounds.csv");
				    dbg << "reaction,LPfile,background_met,volume,baseUB\n";
				    for (auto& r: debugRev)
				        dbg << std::get<0>(r) << ',' << std::get<1>(r) << ',' << std::get<2>(r) << ','
				            << std::get<3>(r) << ',' << std::get<4>(r) << '\n';
				    dbg.close();
				}
		}

		// ──────────────────────────────────────────────────────────────
		// 2) loadAndApplyFBAReactionUpperBounds  (con warning & safe exit su file vuoto)
		// ──────────────────────────────────────────────────────────────
		void loadAndApplyFBAReactionUpperBounds(
				std::vector<FBGLPK::LPprob>& vec_fluxb,
				const std::string& /*unused*/)
		{
				// Trova il primo file esistente nella lista
				const std::string path = firstExisting({
				    "ub_bounds_projected_gui.csv",
				    "ub_bounds_projected.csv",
				    "EX_upper_bounds_FBA.csv"
				});
				if (path.empty()) {
				    std::cerr << "[Warning] FBA bounds CSV not found – UBs unchanged.\n";
				    return;
				}

				// Prova ad aprire il file, esci se non valido
				std::ifstream file(path.c_str());
				if (!file.good()) {
				    std::cerr << "[Warning] Cannot open FBA bounds file '" << path
				              << "' – UBs unchanged.\n";
				    return;
				}

				// Contenitori temporanei
				std::unordered_set<std::pair<std::string,std::size_t>, RxnLPHash> done;
				std::vector<std::tuple<std::string,std::string,double>> debugRows;

				// Salta l'intestazione
				std::string header;
				std::getline(file, header);

				// Leggi e processa ogni riga
				std::string line;
				while (std::getline(file, line)) {
				    if (line.empty()) continue;

				    std::istringstream ss(line);
				    std::string reaction, model, ubStr;
				    std::getline(ss, reaction, ',');
				    std::getline(ss, model,    ',');
				    std::getline(ss, ubStr,     ',');
				    if (reaction.empty() || ubStr.empty()) continue;

				    // Converte upper_bound
				    double newUb;
				    try {
				        newUb = std::stod(ubStr);
				    } catch (...) {
				        std::cerr << "[FBA] Warning: UB non numerico per " << reaction << "\n";
				        continue;
				    }

				    // Trova l'indice LP corrispondente
				    std::size_t lpIdx = findLPIndex(vec_fluxb, model);
				    if (lpIdx == std::numeric_limits<std::size_t>::max()) {
				        std::cerr << "[FBA] Warning: modello '" << model
				                  << "' non trovato per la reazione '" << reaction << "'\n";
				        continue;
				    }

				    auto key = std::make_pair(reaction, lpIdx);
				    if (done.count(key)) continue;
				    done.insert(key);

				    // Recupera la colonna nel problema LP
				    int col = vec_fluxb[lpIdx].fromNametoid(reaction);
				    if (col == -1) {
				        std::cerr << "[FBA] Warning: reazione '" << reaction
				                  << "' non presente in LP '" << vec_fluxb[lpIdx].getFilename()
				                  << "' – skip\n";
				        continue;
				    }

				    // Applica il nuovo upper bound
				    double lb = vec_fluxb[lpIdx].getLwBounds(col);
				    vec_fluxb[lpIdx].update_bound(
				        col,
				        "GLP_DB",
				        trunc(lb, decimalTrunc),
				        trunc(newUb, decimalTrunc)
				    );

				    // Registra per il debug
				    debugRows.emplace_back(reaction, vec_fluxb[lpIdx].getFilename(), newUb);
				}
				file.close();

				// Se nessuna riga valida, esci senza modifiche
				if (debugRows.empty()) {
				    std::cerr << "[Warning] FBA bounds file '" << path
				              << "' is empty or contains no valid entries – UBs unchanged.\n";
				    return;
				}

				// Scrivi il file di debug
				std::ofstream dbg("debug_FBA_updated_bounds.csv");
				dbg << "reaction,LPfile,new_upper_bound\n";
				for (auto& r : debugRows) {
				    dbg << std::get<0>(r) << ','
				        << std::get<1>(r) << ','
				        << std::get<2>(r) << '\n';
				}
				dbg.close();
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
		void updateNonFBAReactionsUB(std::vector<FBGLPK::LPprob>& vec_fluxb,
				                         std::map<std::string,int>&    NumPlaces,
				                         double*                       Value)
		{
				using Row = std::tuple<std::string,std::string,double,double,double,double,double>;
				std::vector<Row> dbg;                                  

				/*──────────────── helper: lista LP contenenti la reazione ───────────────*/
				auto lpListFor = [&](const std::string& rxn)->std::vector<std::size_t>
				{
				    std::vector<std::size_t> v;
				    for (auto& kv : NonFBAReactionBaseUB)
				        if (kv.first.first == rxn) v.push_back(kv.first.second);
				    if (!v.empty()) return v;

				    auto it = reactionToFileMap.find(rxn);
				    if (it != reactionToFileMap.end())
				        v.insert(v.end(), it->second.begin(), it->second.end());
				    return v;
				};

				/*────────── unione delle chiavi di entrambe le mappe (csv + fallback) ───*/
				std::unordered_set<std::string> allRxn;
				for (auto& kv : reactionToFileMap)    allRxn.insert(kv.first);
				for (auto& kv : NonFBAReactionBaseUB) allRxn.insert(kv.first.first);

				/*────────────────────────── loop principale ─────────────────────────────*/
				for (const std::string& rxn : allRxn)
				{
				    /* skip se è una reazione pure-FBA oppure una _f (fissa) */
				    if (FBAreactions.count(rxn) || !endsWith(rxn,"_r"))
				        continue;

				    /* LP che contengono questa reazione */
				    const auto lpIdxList = lpListFor(rxn);
				    if (lpIdxList.empty()) continue;

				    for (std::size_t idx : lpIdxList)
				    {
				        int col = vec_fluxb[idx].fromNametoid(rxn);
				        if (col == -1) continue;

				        /* pop & biomass correnti della specie del LP -------------------*/
				        double pop = 1.0, biomass = 1.0;
				        if (problemBacteriaPlace.count(idx))
				            pop = std::floor(Value[ NumPlaces.at(problemBacteriaPlace[idx]) ]);
				        if (problemBiomassPlace.count(idx))
				            biomass = trunc(Value[ NumPlaces.at(problemBiomassPlace[idx]) ],
				                            decimalTrunc);

				        if (pop <= 0.0) {                              // specie estinta
				            vec_fluxb[idx].update_bound(col,"GLP_FX",0.0,0.0);
				            dbg.emplace_back(rxn,vec_fluxb[idx].getFilename(),
				                             pop,biomass,0.0,
				                             vec_fluxb[idx].getUpBounds(col),0.0);
				            continue;
				        }

				        /* Cₘ (baseUB) ---------------------------------------------------*/
				        auto key = std::make_pair(rxn, idx);
				        auto itB = NonFBAReactionBaseUB.find(key);
				        if (itB == NonFBAReactionBaseUB.end()) continue;
				        const double C_m = fabs(itB->second);          // mmol/L  (o mmol)

				        /* recupero subtype ("exchange" oppure interno) -----------------*/
				        std::string subtype = "exchange";              // fallback
				        auto itSub = ReactionSubtype.find(key);
				        if (itSub != ReactionSubtype.end()) subtype = itSub->second;

				        /* λ -------------------------------------------------------------*/
				        double lambda = 0.0;
				        if (subtype == "exchange") {
				            double denom = 0.0;
				            for (auto k : lpIdxList) {
				                double Nk = std::floor(Value[ NumPlaces.at(problemBacteriaPlace[k]) ]);
				                double Bk = trunc(Value[ NumPlaces.at(problemBiomassPlace[k]) ],
				                                  decimalTrunc);
				                denom += Nk * Bk * 1e-12;              // gDW
				            }
				            lambda = (denom > 0.0) ? 1.0 / denom : 0.0;
				        } else {                                      // internal
				            double denom = pop * biomass * 1e-12;
				            lambda = (denom > 0.0) ? 1.0 / denom : 0.0;
				        }

				        /* nuovo upper-bound --------------------------------------------*/
				        const double newUb = C_m * lambda;            // mmol/gDW/h
				        const double oldUb = vec_fluxb[idx].getUpBounds(col);

				        if (newUb <= 0.0) {
				            vec_fluxb[idx].update_bound(col,"GLP_FX",0.0,0.0);
				        } else {
				            vec_fluxb[idx].update_bound(col,"GLP_DB",0.0,newUb);
				        }

				        dbg.emplace_back(rxn, vec_fluxb[idx].getFilename(),
				                         pop, biomass, C_m, oldUb, newUb);
				    }
				}

				/* dump debug -----------------------------------------------------------*/
				std::ofstream out("debug_nonFBA_runtime_bounds.csv");
				out << "reaction,LPfile,pop,biomass,Cm(oldBase),oldUB,newUB\n";
				for (auto& r : dbg)
				    out << std::get<0>(r) << ','
				        << std::get<1>(r) << ','
				        << std::get<2>(r) << ','
				        << std::get<3>(r) << ','
				        << std::get<4>(r) << ','
				        << std::get<5>(r) << ','
				        << std::get<6>(r) << '\n';
				out.close();
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

				    /* ---------- calcolo nuovo UB secondo Eq.(5.16) ---------- */
					/* ---------- calcolo nuovo UB secondo Eq.(5.16) ---------- */
					double newUB = 0.0;
					double mu = muMaxMap.count(problemIndex) ? muMaxMap[problemIndex] : 1.0;

					if (xB < BioMin) {
							newUB = 0.0;                                  // metabolic distress
							problemsWithLowBiomass[problemIndex] = true;
					}
					else if (Gamma > 0.0) {
							newUB = Gamma / BioMax * mu;                  // mmol / gDW / h
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

};
