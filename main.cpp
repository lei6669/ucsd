#include<vector>
#include<map>
#include<iostream>
#include<iterator>
#include<algorithm>
#include<chrono>
#include<fstream>
#include <sstream>
#include <omp.h>
#include <unistd.h>

int THREAD_NUM = 8;

typedef unsigned int SID;
typedef unsigned int Intensity;
typedef unsigned int MZ;

typedef std::pair<MZ, Intensity> Peak;
typedef std::vector<Peak> Spectrum;
typedef std::vector<Spectrum> RawData;

// typedef std::pair<SID, Intensity> BucketPeak;
typedef std::pair<SID, std::pair<MZ, Intensity> > BucketPeak;

typedef std::vector<BucketPeak> Bucket;

typedef std::vector<Bucket> Index;

// static const MZ MAX_MZ = 20000;
static const MZ MAX_MZ = 1000;
// static const int num_buckets = 1; //# of bins per mz used for index
static const int num_buckets = 20;

RawData * load_raw_data(char *file, int &total_spectra, int &num_peaks) {
	RawData * spectra = new RawData();
	unsigned int file_id;
	unsigned int num_spectra;           //total number of spectra inside the file
	std::vector< std :: pair<unsigned int, unsigned int> > position; // starting position and ending position for each spectrum
	std::ifstream in(file, std::ios::in | std::ios::binary);
	in.read((char*)&file_id, sizeof( unsigned int ));
	in.read((char*)&num_spectra, sizeof( unsigned int ));
	position.resize(num_spectra);

	//total_spectra = 0 as input means to load all spectra from the DB
	if (total_spectra == 0){
		total_spectra = num_spectra;
	}

	if (num_spectra > 0) {
		position[0].first = num_spectra + 2; //starting offset in the file of the first spectrum
		in.read((char *) &position[0].second, sizeof(unsigned int)); //ending position of the first spectrum

		for (unsigned int spec_idx = 1; spec_idx < num_spectra; ++spec_idx) {
			position[spec_idx].first = position[spec_idx - 1].second;  //starting position
			in.read((char *) &position[spec_idx].second, sizeof(unsigned int)); //ending position
		}

		for (unsigned int spec_idx = 0; spec_idx < total_spectra; ++spec_idx) {
			in.seekg(position[spec_idx].first * 4 + 16); //setting the offset to read the spectrum
			unsigned int size = position[spec_idx].second - position[spec_idx].first - 4; //total number of peaks inside the spectrum

			//num_peaks = 0 as input means to load all peaks of the spectrum
			if (num_peaks == 0 || num_peaks < size){
				num_peaks = size;
			}

			Spectrum spectrum;
			//Populating peak info per spectrum
			for (int i = 0; i < num_peaks; ++i) {
				unsigned int peak;
				unsigned int mz;
				in.read((char *) &peak, sizeof(unsigned int));
				mz = peak >> 8;
				spectrum.push_back(Peak(mz,peak - (mz<<8)));
			}
			spectra -> push_back(spectrum);
		}
	}
	in.close();
	return spectra;
}

void dump_spectrum(Spectrum *s) {
	std::cerr << "[";
	for(auto & p: *s) {
		std::cerr << "[" << p.first << ", " << p.second << "],";
	}
	std::cerr << "]";
}

// *helper fuction for printing
void dump_bucket(Bucket *bucket) {
	std::cerr << "[";
	for(auto & p: *bucket) {
		std::cerr << "[" << p.first << ", " << p.second.second << "],";
	}
	std::cerr << "]";
}

void dump_raw_data(RawData* r) {

	int c = 0;
	for(auto & s: *r) {
		std::cerr << "SID=" << c;
		dump_spectrum(&s);
		c++;
		std::cerr << "\n";
	}
}

void dump_index(Index *index) {
	for(MZ mz = 0; mz < MAX_MZ; mz++) {
		if(!(*index)[mz].empty()) {
			std::cerr << "MZ = " << mz << ": ";
			dump_bucket(&(*index)[mz]);
			// dump_spectrum(&(*index)[mz]);
			std::cerr << "\n";
		}
	}
}

void json_reconstruction(char * file, const std::vector<std::map<SID, Spectrum>> &reconstructed_spectra) {

	std::ofstream out(file, std::ios::out);


	//auto end = reconstructed_spectra.end();
	out << "[\n";
	for(const auto & queries_results: reconstructed_spectra) {
		out << "[\n";
		for(const auto &[sid, spectrum]:queries_results) {

			out << "\t{ " << "";
			out << "\"" << sid << "\": " << "[\n";
			for (auto j = spectrum.begin(); j != spectrum.end(); j++) {
				out << "\t\t\t[ " << j->first << ", " << j->second << " ]";
				if (std::next(j) != spectrum.end()) {
					out<<  ",";
				}
				out << "\n";
			}
			out << "\t\t]\n";
			out << "\t}";

			if (&sid != &queries_results.rbegin()->first) {
				out << ",";
			}

			out << "\n";
		}
		out << "]";
		if (&queries_results != &*reconstructed_spectra.rbegin()) {
			out << ",";
		}
		out << "\n";
	}
	out << "]\n";
}

Spectrum * load_query(char*file) {
	Spectrum *n = new Spectrum;
	std::ifstream in(file);
	std::string line;

	while (std::getline(in, line)) {
		std::vector<unsigned int> lineData;
		std::stringstream lineStream(line);
		unsigned int value;
		while (lineStream >> value) {
			lineData.push_back(value);
		}
		n->push_back(Peak(lineData[0], lineData[1]));
	}
	return n;
}

std::vector<Spectrum> * load_queries(char*file) {
	auto n = new std::vector<Spectrum>;

	std::ifstream in(file);
	std::string line;

	int spectra_count = 0;
	in >> spectra_count;
	for(int i = 0; i < spectra_count; i++) {
		Spectrum s;
		int peak_count;
		in >> peak_count;
		for (int p = 0; p < peak_count; p++) {
			MZ mz;
			Intensity intensity;
			in >> mz >> intensity;
			s.push_back(Peak(mz, intensity));
		}
		n->push_back(s);
	}
	return n;
}

Index * build_index(RawData * data) {
	// set up number of threads
	omp_set_num_threads(THREAD_NUM);
	Index *index = new Index;
	// initiate all buckets
	for(MZ mz = 0; mz < MAX_MZ; mz++) {
		index->push_back(Bucket());
	}
	// create mutex lock for each bucket
	omp_lock_t index_lock[MAX_MZ];
	for (int i=0; i<MAX_MZ; i++) {
		omp_init_lock(&index_lock[i]);
	}
	auto reconstruct_start = std::chrono::high_resolution_clock::now();
	// multi-thread for building index on each raw data spectrum
	int size = data->size();
	//#pragma omp parallel for
	for(SID sid = 0; sid < size; sid++) {
		unsigned int unit_frag;
		Spectrum spectrum = (*data)[sid];
		// for(auto & peak: (*data)[sid]) {
		for(int i=0; i<spectrum.size(); i++) {
			Peak peak = spectrum[i];
			unit_frag = peak.first/num_buckets;
			if (unit_frag < MAX_MZ) {
				// (*index)[unit_frag].push_back(BucketPeak(sid, std::pair<MZ, Intensity>(peak.first, peak.second)));
				// (*index)[unit_frag].push_back(BucketPeak(sid, peak.second));
				// omp_set_lock(&index_lock[unit_frag]);
				int currentSize = (*index)[unit_frag].size();
				if (currentSize == 0 || (*index)[unit_frag][currentSize-1].first < sid) {
					// std::cout<<"aa"<<std::endl;
					//omp_set_lock(&index_lock[unit_frag]);
					(*index)[unit_frag].push_back(BucketPeak(sid, std::pair<MZ, Intensity>(peak.first, peak.second)));
					//omp_unset_lock(&index_lock[unit_frag]);
				} else {  // now need to push the element into the front of equal sid
					// std::cout<<"ss"<<std::endl;
					if ((*index)[unit_frag][currentSize-1].second.second <= peak.second) {
						//omp_set_lock(&index_lock[unit_frag]);
						(*index)[unit_frag].push_back(BucketPeak(sid, std::pair<MZ, Intensity>(peak.first, peak.second)));
						//omp_unset_lock(&index_lock[unit_frag]);
					} else {
						int insertIndex = currentSize-1;
						while (insertIndex-1>=0 && (*index)[unit_frag][insertIndex-1].first == sid && (*index)[unit_frag][insertIndex-1].second.second>peak.second) {
							insertIndex -= 1;
						}
						//omp_set_lock(&index_lock[unit_frag]);
						(*index)[unit_frag].insert((*index)[unit_frag].begin()+insertIndex, BucketPeak(sid, std::pair<MZ, Intensity>(peak.first, peak.second)));
						//omp_unset_lock(&index_lock[unit_frag]);
					}
				}
				// omp_unset_lock(&index_lock[unit_frag]);
			}
		}
	}
	//#pragma omp barrier
	auto reconstruct_end = std::chrono::high_resolution_clock::now();
	// multi-thread on sroting part
	// #pragma omp parallel for
	// for(MZ mz = 0; mz < MAX_MZ; mz++) {
	// 	std::sort((*index)[mz].begin(), (*index)[mz].end());
	// }

	// int who = 999;
	// int compare[(*index)[who].size()];
	// int inten[(*index)[who].size()];
	// for (int i=0; i<(*index)[who].size(); i++) {
	// 	compare[i] = (*index)[who][i].first;
	// 	inten[i] = (*index)[who][i].second.second;
	// }
	// std::sort((*index)[who].begin(), (*index)[who].end());
	// #pragma omp barrier
	auto reconstruct_ended = std::chrono::high_resolution_clock::now();
	std::cerr << "     1:" << (std::chrono::duration_cast<std::chrono::nanoseconds>(reconstruct_end - reconstruct_start).count()+0.0)/1e9 << " s\n";
	std::cerr << "     2:" << (std::chrono::duration_cast<std::chrono::nanoseconds>(reconstruct_ended - reconstruct_end).count()+0.0)/1e9 << " s\n";
	// std::cout<< (*index)[999][0].first<<" "<< (*index)[999][1].first<<" "<< (*index)[999][2].first<<" "<< (*index)[999][3].first<<" "<< (*index)[999][4].first <<std::endl;

	// for (int i=0; i<(*index)[who].size(); i++) {
	// 	if ((*index)[who][i].first != compare[i] || (*index)[who][i].second.second != inten[i]) {
	// 		std::cout<< (*index)[who][i].first <<" "<< (*index)[who][i].second.second<< " "<< compare[i] << " "<< inten[i]<< std::endl;
	// 	}
	// 	// std::cout<< (*index)[who][i].first << " "<< compare[i] << std::endl;
	// }
	// exit(2);
	return index;
}

std::vector<std::map<SID, Spectrum>> *reconstruct_candidates(Index * index, const std::vector<Spectrum> & queries) {
	omp_lock_t map_lock;
	omp_init_lock(&map_lock);
	auto reconstructed_spectra = new std::vector<std::map<SID, Spectrum>>;
	std::map<int, std::map<SID, Spectrum>> temp;
	int size = queries.size();
	#pragma omp parallel for
	for(int i=0; i<size; i++) {
		Spectrum query = queries[i];
		std::map<SID, Spectrum> m;
		for(auto & query_peak: query) {
			unsigned int unit_q_mz = query_peak.first / num_buckets;
			if ((*index).size() > unit_q_mz) {
				for(auto & bucket_peak : (*index)[unit_q_mz]) {
					if (bucket_peak.second.first == query_peak.first) {
						// m[bucket_peak.first].push_back(Peak(query_peak.first, bucket_peak.second));
						m[bucket_peak.first].push_back(Peak(query_peak.first, bucket_peak.second.second));
					}
				}
			}
		}
		omp_set_lock(&map_lock);
		temp[i] = m;
		omp_unset_lock(&map_lock);
	}
	for (int i=0; i<size; i++) {
		reconstructed_spectra->push_back(temp[i]);
	}
	return reconstructed_spectra;
}

int main(int argc, char * argv[]) {

	if (argc != 4) {
		std::cerr << "Usage: main <raw data file> <query file> <output json>\n";
		exit(1);
	}

	/* Declare number of spectra you want to load from the database
	 * 0 - load all spectra from the file
	 * n - load first N spectra from the file
	 * for initial test in example n = 3
	 */
	int total_spectra = 0;

	/* Declare number of peaks you want to load from each spectrum
	 * 0 - load all peaks from the spectrum
	 * m - load first N peaks from the spectrum
	 * for initial test in example m = 5
	 */
	int num_peaks = 0;

	bool demo = false;

	if (std::getenv("DEMO")) {
		total_spectra = 3;
		num_peaks = 5;
		demo = true;
	}

	RawData * raw_data = load_raw_data(argv[1], total_spectra, num_peaks);
	//std::cerr << "raw_data=\n";
	//dump_raw_data(raw_data);

	std::vector<Spectrum> *queries = load_queries(argv[2]);
	if (demo) {
		std::cerr << "queries=\n";
		for(auto &s: *queries) {
			dump_spectrum(&s);
			std::cerr << "\n";
		}
		std::cerr << "\n";
	}

	// Here's where the interesting part starts

	auto index_build_start = std::chrono::high_resolution_clock::now();
	Index * index = build_index(raw_data);
	auto index_build_end = std::chrono::high_resolution_clock::now();

	delete raw_data;

	if (demo) {
		std::cerr << "Index: ";
		dump_index(index);
	}

	auto reconstruct_start = std::chrono::high_resolution_clock::now();
	auto reconstructed_spectra = reconstruct_candidates(index, *queries);
	auto reconstruct_end = std::chrono::high_resolution_clock::now();
	json_reconstruction(argv[3], *reconstructed_spectra);

	std::cerr << "Found " << reconstructed_spectra->size() << " candidates \n";
	std::cerr << "Building the index took " << (std::chrono::duration_cast<std::chrono::nanoseconds>(index_build_end - index_build_start).count()+0.0)/1e9 << " s\n";
	std::cerr << "Reconstruction took     " << (std::chrono::duration_cast<std::chrono::nanoseconds>(reconstruct_end - reconstruct_start).count()+0.0)/1e9 << " s\n";
}
