#include "General.h"

int ImagingAlgorithm::ErrorInt()
{
    static int fErrorInt = -999999;
    return fErrorInt;
}

double ImagingAlgorithm::ErrorDouble()
{
    static double fErrorDouble = -1e10;
    return fErrorDouble;
}

TChain *ImagingAlgorithm::GenerateChain(std::string treeName, const std::vector<std::string> &keyWord, TChain *chain, string sFolder, vector<std::string> notInclude)
{
    gSystem->Exec(Form("ls %s/ 1> .~filelist 2> /dev/null", sFolder.c_str()));

    auto ch = chain;
    if (ch == NULL)
    {
        ch = new TChain(treeName.c_str());
    }

    ifstream file_list(".~filelist");
    for (int i = 0; file_list.is_open() && file_list.eof() == false; i++)
    {
        string s_temp;
        file_list >> s_temp;

        bool continueFlag = 0;
        for (size_t keyCount = 0; keyCount < keyWord.size(); keyCount++)
        {
            if (s_temp.find(keyWord[keyCount]) == string::npos)
            {
                continueFlag = 1;
            }
            else
            {
                // find in notInclude vector
                for (size_t i = 0; i < notInclude.size(); i++)
                {
                    if (s_temp.find(notInclude[i]) != string::npos)
                    {
                        continueFlag = 1;
                    }
                }
            }
        }
        if (continueFlag)
            continue;

        cout << "File: " << s_temp << " Read" << endl;
        string fileName = sFolder + "/" + s_temp;
        // ch->Add(s_temp.c_str());
        ch->Add(fileName.c_str());
    }
    cout << "Totally get " << ch->GetEntries() << " Entries" << endl;
    gSystem->Exec("rm .~filelist");
    return ch;
}
