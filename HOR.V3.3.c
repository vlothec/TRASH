#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// TODO: change alilength calculation to let fasta alignments that are wrapped to be used as input (line 65)
char **allAlignedSeqs;
int alilength = 0;
int threshold;          //less OR equal than threshold variants allowed
int cutoff;             //cutoff is the smallest detectable HOR (smaller are removed, equal are kept)
int split;

bool compareAB(int *indA, int *indB, int *snvCounto);


/*
args:
0 exe
1 location
2 alignment file name
3 split
4 threshold
5 cutoff
6 method
*/


int main(int argc, char *argv[])
{
    printf( "***********************************************************************\n"
            "*                     HOR identification software                     *\n"
            "*                                V3.3                                 *\n"
            "*        Single file input split into repeats from two regions        *\n"
            "*  Repeats need to have their direction set in the name as D1 or D2!  *\n"
            "***********************************************************************\n\n\n");
    clock_t t;
    t = clock();

    if(argc == 7)
    {
        printf("Comparing repeats from the folder: %s\nfile: %s, repeats split after repeat no: %s, threshold: %s and cutoff value: %s type is %s \n\n", argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
    }
    else
    {
        printf("Wrong amount of arguments supplied.\nProvide (in order):\nDirectory of .fasta alignments to be compared and output placed\nName of the alignment\nWhere the alignment should be split\nThreshold value\nCutoff value\nComparison type: 1(individual) or 2(double)\n\n");
        return 1;
    }


    split = atoi(argv[3]) - 1;
    threshold = atoi(argv[4]);
    cutoff = atoi(argv[5]);

    FILE *log;                                                                                                          // Save the output
    char logfile[1000];
    sprintf(logfile,"%slog_%s_t_%i_c_%i.txt", argv[1], argv[2], threshold, cutoff);
	log = fopen(logfile, "w");

    char alifile[1000];

    strcpy(alifile, argv[1]);
    strcat(alifile, argv[2]);


    printf("Alifile: %s\n", alifile);
    // Count how many sequences and what is the length of the alignment
    printf("Checking the alignment file and calculating the number of sequences ...\n\n");

    FILE *alignment;
	alignment = fopen(alifile, "r");

	if(alignment == NULL)
    {
        printf("Error: the file does not exist, is corrupted, or denied access\n\n");
        exit(0);
    }

    for(int c = fgetc(alignment); c != 10; c = fgetc(alignment)) {   }
    for(int c = fgetc(alignment); c != '>'; c = fgetc(alignment))
    {
        fprintf(log, "%c ", c);
        if(c != 10)
        {
            alilength = alilength + 1;
        }
    }
    printf("The length of the alignment is: %i\n\n", alilength);



	int alicount = 2;
	for (int c = fgetc(alignment); c != EOF; c = fgetc(alignment))
	{
		if(c == '>')
		{
			alicount = alicount + 1;
		}
	}
	printf("Count of aligned sequences is: %i\n\n", alicount);

	fclose(alignment);



    alignment = fopen(alifile, "r");                                                                                // read aligned sequences into a new allAlignedSeqs string array

    allAlignedSeqs = malloc(alicount * sizeof(char*));
    printf("Allocating memory for the alignment ...\n\n");
    int dir[alicount];
    int i,j,k;
    int readPos = 0;
    char c;

    char temp[alilength+2];
    for(i = 0; i < alicount; i++)
    {
        allAlignedSeqs[i] = malloc((alilength + 2) * sizeof(char));
    }
    printf("Reading the alignment ...\n\n");
	for(j = 0; j < alicount; j ++)
	{
        for(k=0; k < alilength; k++)
        {
            temp[k] = 0;
        }
		readPos = 0;
        for(c = fgetc(alignment); c != 'D'; c = fgetc(alignment))
		{
		    //printf("%c", c);
        }
        dir[j] = 0;
        dir[j] = fgetc(alignment);
		for(c = fgetc(alignment); c != 10; c = fgetc(alignment))
		{

        }
		for(c = fgetc(alignment); ((c != '>') && (c != EOF)); c = fgetc(alignment))
		{
		    if(c != 10)
            {
                temp[readPos] = c;
                readPos++;
            }

		    //printf("%c", temp[readPos]);
			//readPos++;
		}
		readPos++;
		temp[readPos] = '\0';
		//printf("\n%s", temp);
        //printf("len of allAlignedSeqs[0] is %i\n", strlen(allAlignedSeqs[0]));

        strcpy(allAlignedSeqs[j], temp);
        //printf("%s", temp);

	}
    fclose(alignment);
    printf("closing alignment \n\n");

    for(int i = 0; i < alilength; i++)
    {
        //printf("2dir %c\n", dir[i]);
    }



    printf("\nprepare the output file...\n\n");
    //prepare the output file
    FILE *csvfile;                                                                                                          // Save the output
    char csvfilename[1000];
    sprintf(csvfilename,"%sHORs_method_%i_%s_t_%i_c_%i.csv", argv[1], atoi(argv[6]), argv[2], threshold, cutoff);
	csvfile = fopen(csvfilename, "w");
    fprintf(csvfile, "start_A,end_A,start_B,end_B,direction(1=para_2=perp),total_variant\n");
	//later do:
	//fprintf(csvfile, "%i,%i,%i,%i,%i,%i,%i\n", i+1, horSTARTAdata, horENDAdata, horSTARTBdata, horENDBdata, horDIRdata,horSNVdata);

	//prepare variables to store HOR stats currently processed

    int snvcount;
    int openHOR = 0;
    int snvBack = 0;
    bool isSimilar = 0;

    fprintf(log, "Starting main analysis...\n\n");
    if(atoi(argv[6]) == 1)
    {
        fprintf(log, "\nAnalysing all repeats as a whole (without dividing)\n");

        fprintf(log, "\nhalf I(like 2.III): ");
        for(int i = 0; i < alicount/100; i++)
        {
            fprintf(log, "o");
        }
        fprintf(log, "o\nhalf I(like 2.III): ");

        for(int s = alicount - 1; s > 0; s--)
        {
            snvcount = 0;
            openHOR = 0;
            i = 0;
            j = s;
            if(j%100 == 0)
            {
                fprintf(log, "I");
            }
            while((j < alicount) && (i < alicount))
            {
                //printf("i %i, j %i\n", i, j);
                snvBack = snvcount;
                isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] == dir[j]);

                //printf (" isSimilar is %i, i is %i, j is  %i, snvcount is %i, openHOR is  %i", isSimilar, i, j, snvcount, openHOR);
                if(isSimilar && (openHOR == 0))
                {
                    //printf(" open loop");
                    openHOR = 1;
                } else if(isSimilar && (openHOR > 0))
                {
                    openHOR = openHOR + 1;
                } else if(!isSimilar && (openHOR > 0))
                {
                    //printf(", close loop");
                    if(openHOR >= cutoff)
                    {
                        //printf(", and save HOR: %i,%i,%i,%i,%i,%i", i - openHOR + 1, i, j - openHOR + 1, j, 1, snvBack);
                        fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j - openHOR + 1, j, 1, snvBack);
                    }
                    openHOR = 0;
                    snvcount = 0;
                }
                i++;
                j++;
            }
            if(isSimilar && (openHOR >= cutoff))
            {
                //printf(", and save HOR at the end: %i,%i,%i,%i,%i,%i", i - openHOR + 1, i, j - openHOR + 1, j, 1, snvBack);
                fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j - openHOR + 1, j, 1, snvBack);
            }
            isSimilar = 0;
        }




        fprintf(log, "\n\nquadrant II(like 2.I): ");

        for(int i = 0; i < alicount/100; i++)
        {
            fprintf(log, "o");
        }
        fprintf(log, "o\nquadrant II(like 2.I): ");

        for(int s = 1; s < (alicount - 1) ; s++)
        {
            snvcount = 0;
            openHOR = 0;
            i = 0;
            j = s;
            if(j%100 == 0)
            {
                fprintf(log, "I");
            }
            while(j > i)
            {
                //printf("i %i, j %i\n", i, j);
                snvBack = snvcount;
                isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] != dir[j]);

                //printf (" isSimilar is %i, i is %i, j is  %i, snvcount is %i, openHOR is  %i", isSimilar, i, j, snvcount, openHOR);
                if(isSimilar && (openHOR == 0))
                {
                    //printf(", open loop");
                    openHOR = 1;
                } else if(isSimilar && (openHOR > 0))
                {
                    //printf(" continue HOR");
                    openHOR = openHOR + 1;
                } else if(!isSimilar && (openHOR > 0))
                {
                    //printf(", close loop");
                    if(openHOR >= cutoff)
                    {
                        //printf(", and save HOR: %i,%i,%i,%i,%i,%i", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
                        fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
                    }
                    openHOR = 0;
                    snvcount = 0;
                }
                i++;
                j--;
            }
            //printf(", openHOR is %i and isSimilar is %i", openHOR, isSimilar);
            if(isSimilar && (openHOR >= cutoff)) // close an open HOR if it was opened but not closed due to running into the edge of the array
            {
                //printf(", and save HOR at the end: %i,%i,%i,%i,%i,%i", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
                fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
            }
            isSimilar = 0;
        }


        snvcount = 0;
        openHOR = 0;
        fprintf(log, "\nquadrant III(like 2.II): ");
        for(int i = 0; i < alicount/100; i++)
        {
            fprintf(log, "o");
        }
        fprintf(log, "\nquadrant III(like 2.II): ");
        for(int s = 0; s < (alicount - 1); s++)
        {
            snvcount = 0;
            openHOR = 0;
            i = s;
            j = alicount - 1;
            if(i%100 == 0)
            {
                fprintf(log, "I");
            }

            while(j > i)
            {
                //printf("i %i, j %i\n", i, j);
                snvBack = snvcount;
                isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] != dir[j]);

                //printf (" isSimilar is %i, i is %i, j is  %i, snvcount is %i, openHOR is  %i", isSimilar, i, j, snvcount, openHOR);
                if(isSimilar && (openHOR == 0))
                {
                    //printf("open loop");
                    openHOR = 1;
                } else if(isSimilar && (openHOR > 0))
                {
                    //printf(" continue HOR");
                    openHOR = openHOR + 1;
                } else if(!isSimilar && (openHOR > 0))
                {
                    //printf(", close loop");
                    if(openHOR >= cutoff)
                    {
                        //printf(", and save HOR: %i,%i,%i,%i,%i,%i", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
                        fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
                    }
                    openHOR = 0;
                    snvcount = 0;
                }
                i++;
                j--;
            }
            //printf(", openHOR is %i and isSimilar is %i", openHOR, isSimilar);
            if(isSimilar && (openHOR >= cutoff)) // close an open HOR if it was opened but not closed due to running into the edge of the array
            {
                //printf(", and save the end line HOR: %i,%i,%i,%i,%i,%i", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
                fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
            }
            isSimilar = 0;
        }



        fprintf(log, "\n");



    } else
    {
        fprintf(log, "\n\nquadrant I: ");


        for(int i = 0; i < (alicount - split)/100; i++)
        {
            fprintf(log, "o");
        }
        fprintf(log, "o\nquadrant I: ");

        for(int s = split+1; s < alicount ; s++)
        {
            snvcount = 0;
            openHOR = 0;
            i = 0;
            j = s;
            if(j%100 == 0)
            {
                fprintf(log, "I");
            }
            while((j > split) && (i <= split))
            {
                //printf("i %i, j %i\n", i, j);
                snvBack = snvcount;
                isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] != dir[j]);

                if(isSimilar && (openHOR == 0))
                {
                    //printf("open loop");
                    openHOR = 1;
                } else if(isSimilar && (openHOR > 0))
                {
                    openHOR = openHOR + 1;
                } else if(!isSimilar && (openHOR > 0))
                {
                    if(openHOR >= cutoff)
                    {
                        fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
                    }
                    openHOR = 0;
                    snvcount = 0;
                }
                i++;
                j--;
            }
            if(isSimilar && (openHOR >= cutoff)) // close an open HOR if it was opened but not closed due to running into the edge of the array
            {
                fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
            }
            isSimilar = 0;
        }


        fprintf(log, "\nquadrant II: ");
        for(int i = 0; i < split/100; i++)
        {
            fprintf(log, "o");
        }
        fprintf(log, "\nquadrant II: ");
        for(int s = 1; s <= split; s++)
        {
            snvcount = 0;
            openHOR = 0;
            i = s;
            j = alicount - 1;
            if(i%100 == 0)
            {
                fprintf(log, "I");
            }

            while((j > split) && (i <= split))
            {
                //printf("i %i, j %i\n", i, j);
                snvBack = snvcount;
                isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] != dir[j]);


                if(isSimilar && openHOR == 0)
                {
                    //printf("open loop");
                    openHOR = 1;
                } else if(isSimilar && (openHOR > 0))
                {
                    openHOR = openHOR + 1;
                } else if(!isSimilar && (openHOR > 0))
                {
                    if(openHOR >= cutoff)
                    {
                        fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
                    }
                    openHOR = 0;
                    snvcount = 0;
                }
                i++;
                j--;
            }
            if(isSimilar && (openHOR >= cutoff)) // close an open HOR if it was opened but not closed due to running into the edge of the array
            {
                fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j + 2, j + openHOR + 1, 2, snvBack);
            }
            isSimilar = 0;
        }

        fprintf(log, "\nquadrant III: ");
        for(int i = 0; i < (alicount - split)/100; i++)
        {
            fprintf(log, "o");
        }
        fprintf(log, "o\nquadrant III: ");


        for(int s = alicount - 1; s > split; s--)
        {
            snvcount = 0;
            openHOR = 0;
            i = 0;
            j = s;
            if(j%100 == 0)
            {
                fprintf(log, "I");
            }
            while((j < alicount) && (i <= split))
            {
                //printf("i %i, j %i\n", i, j);
                snvBack = snvcount;
                isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] == dir[j]);

                if(isSimilar && (openHOR == 0))
                {
                    //printf("open loop");
                    openHOR = 1;
                } else if(isSimilar && (openHOR > 0))
                {
                    openHOR = openHOR + 1;
                } else if(!isSimilar && (openHOR > 0))
                {
                    if(openHOR >= cutoff)
                    {
                        fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j - openHOR + 1, j, 1, snvBack);
                    }
                    openHOR = 0;
                    snvcount = 0;
                }
                i++;
                j++;
            }
            if(isSimilar && (openHOR >= cutoff))
            {
                fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j - openHOR + 1, j, 1, snvBack);
            }
            isSimilar = 0;
        }

        fprintf(log, "\nquadrant IV: ");
        for(int i = 0; i < split/100; i++)
        {
            fprintf(log, "o");
        }
        fprintf(log, "\nquadrant IV: ");
        for(int s = 1; s <= split; s++)
        {
            snvcount = 0;
            openHOR = 0;
            i = s;
            j = split + 1;
             if(i%100 == 0)
            {
                fprintf(log, "I");
            }
            while((j < alicount) && (i <= split))
            {
                //printf("i %i, j %i\n", i, j);
                snvBack = snvcount;
                isSimilar = compareAB(&i, &j, &snvcount) && (dir[i] == dir[j]);

                if(isSimilar && (openHOR == 0))
                {
                    //printf("open loop");
                    openHOR = 1;
                } else if(isSimilar && (openHOR > 0))
                {
                    openHOR = openHOR + 1;
                } else if(!isSimilar && (openHOR > 0))
                {
                    if(openHOR >= cutoff)
                    {
                        fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j - openHOR + 1, j, 1, snvBack);
                    }
                    openHOR = 0;
                    snvcount = 0;
                }
                i++;
                j++;
            }
            if(isSimilar && (openHOR >= cutoff))
            {
                fprintf(csvfile, "%i,%i,%i,%i,%i,%i\n", i - openHOR + 1, i, j - openHOR + 1, j, 1, snvBack);
            }
            isSimilar = 0;
        }
        fprintf(log, "\n");
    }

    fprintf(log, "Finished the main analysis...\n\n");

   for(int i = 0; i < alicount; i++)
    {
        free(allAlignedSeqs[i]);
    }
    fprintf(log, "Freed memory\n");
    free(allAlignedSeqs);

    //printf("Freed the seq memory...\n\n");

    fprintf(log, "\n\nHOR output file created:\n\n%sHORs_%s_t_%i_c_%i.csv\n\n", argv[1], argv[2], threshold, cutoff);
    fclose(csvfile);
    t = clock() - t;


    fprintf(log, "HOR output file created:\n\n%sHORs_%s_t_%i_c_%i.csv\n\n", argv[1], argv[2], threshold, cutoff);
    fprintf(log, "The program finished in %f minutes.\n",((float)t)/CLOCKS_PER_SEC/60);
    fclose(log);
    printf ("The program finished in %f minutes.\n",((float)t)/CLOCKS_PER_SEC/60);

	return 0;

}





bool compareAB(int *indA, int *indB, int *snvCounto)
{
    int snv = 0;
    for(int i = 0; i < alilength; i++)
    {
        snv = snv + (allAlignedSeqs[*indA][i] != allAlignedSeqs[*indB][i]);
        //printf("%c", allAlignedSeqs[*indA][i]);
    }
    //printf("\n");
    if(snv <= threshold)
    {
        *snvCounto = *snvCounto + snv;
    }
    //printf ("\n comparing %i with %i,  snc count is %i", *indA, *indB, snv);
    return(snv <= threshold);
}
