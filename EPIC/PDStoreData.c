/*@Start***********************************************************
 * GEMSBG Source File
 * Copyright (C) 1989-1998 The General Electric Company
 *
 *      File Name:  PDStoreData.c
 *      Developer:  B. L. Mazin
 *
 * $Source: PDStoreData.c $
 * $Revision: 1.5 $  $Date: 6/23/98  19:04:29  $
 *
 *@Description
 *      Store Pulse Data into a file.
 *      This routine reads global memory where the waveform and
 *      and instruction memory is saved, and stores it to disk.
 *
 *
 * do not edit anything above this line
 *
 ******************************************************************
 *      Revision History
 ******************************************************************
 Version          Date/        Author
                  Comment
 ------------------------------------------------------------------
                  05/16/95     JDM
                  Commented out #include <widec.h>, as it multiply defined 
                  some of stdio.h contents.

 sccs1.5          22-Jun-98    Dale Thayer
                  Original for cardiac, CV1, conversion to ANSI C.

 /main/mr_main/3  18-Nov-98    Dale Thayer
                  Original ClearCase consolidated version for ANSI C.

		  03-Mar-01    Angus Bond
		  Convert to PD/PDX File Format, Version 2.0.  For a 
		  description of the file format, see "WTools PD/PDX File 
		  Format", dated 02 March 2001, by Angus Bond.

 *@End*************************************************************/

#include <stdio.h>		/* Defines NULL, ... */
#include <string.h>		/* Defines strlen(), ... */
#include <time.h>		/* Defines localtime(), strftime(), ... */
#include <sys/types.h>		/* Defines time_t, ... */
#include "XmBuilderUtils.h"

#include "stddef_ep.h"
#include "pulsegen.h"		/* Defines WF_MAX_PROCESSORS, ... */
#include "ppulsedef.h"
#include "plotdata.h"		/* Defines PDPlotDataHeader, ... */
#if defined(IPG_TGT)
#include "pgen_rsp_glob.h"
#elif defined(MGD_TGT)
#include "pgen_globals.h"	/* Defines NUM_WAVEGEN, NUM_DSP, rsp_p, ... */
#if (NUM_WAVEGEN != WF_MAX_PROCESSORS)
#error NUM_WAVEGEN and WF_MAX_PROCESSORS do not have the same value.
#endif
#endif
#include "psdplotlib_proto.h"
#include "psdutil.h"		/* Defines OpenFile(), ... */
#include "epic.global.h"
#include "epicfuns.h"

/*****************************************************************************/
/*
 * Local Function Prototypes...
 */

static CHAR *constructPDFileName( CHAR *rootfile );
static CHAR *constructPSDName( CHAR *rootfile );

/*****************************************************************************/
/*
 * This routine will store the data from the pulsegen host file to
 *  a database
 */

STATUS 
PDStoreData( CHAR *rootfile )
{
    CHAR *filename;
    CHAR *psdname;
    FILE *plotdata_ptr;			/* File Ptr */
    const CHAR *cdata;
    CHAR cbuf[80];
    time_t t;
    SHORT s;
    CHAR *cp;
    INT i;
    LONG jj;
#if defined(IPG_TGT)
    PDIpgPlotDataHeader plotdata_hdr;	/* Hdr Packet */
    LONG base;
    LONG end;
    LONG bytes;
    SHORT *wmem;
    SHORT *imem;
#elif defined(MGD_TGT)
    n32 n32_datum;
    void *pointer_datum;
    SeqInstruction0_s *imem;
    n16 *wmem;
#endif

    /**************************************/
    /******** Open the PD/PDX file ********/
    /**************************************/

    if ((filename = constructPDFileName(rootfile)) == NULL)
        return FAILURE;
    if ((psdname = constructPSDName(rootfile)) == NULL)
        return FAILURE;
#if (WTOOLS_DIAG >= 1)
    printf("PDStoreData(): rootfile = \"%s\"\n", rootfile);
    printf("PDStoreData(): filename = \"%s\"\n", filename);
    printf("PDStoreData(): psdname  = \"%s\"\n", psdname);
    fflush(stdout);
#endif

    if ((plotdata_ptr = OpenFile(filename, WRITE_MODE)) == NULL)
    {
        WTFree(filename);
        WTFree(psdname);
        return FAILURE;
    }

    /**************************************/
    /****** Write the header section ******/
    /**************************************/

    cdata = "PD/PDX File\n";
    fwrite(cdata, sizeof(CHAR), (strlen(cdata) + 1), plotdata_ptr);

#if defined(IPG_TGT)
    cdata = "IPG Legacy\n";
#elif defined(MGD_TGT)
    cdata = "LxMGD 1.0\n";
#endif
    fwrite(cdata, sizeof(CHAR), (strlen(cdata) + 1), plotdata_ptr);

    t = time(NULL);
    strftime(cbuf, (sizeof(cbuf) - 1), "Created %Y-%m-%d %H:%M:%S\n", 
             localtime(&t));
    fwrite(cbuf, sizeof(CHAR), (strlen(cbuf) + 1), plotdata_ptr);

    s = 0x1122;
    cp = (CHAR *) &s;
    if (*cp == 0x11)
        cdata = "Big-endian\n";
    else if (*cp == 0x22)
        cdata = "Little-endian\n";
    else
        cdata = "Indeterminate endian-ness\n";
    fwrite(cdata, sizeof(CHAR), (strlen(cdata) + 1), plotdata_ptr);

    cdata = "End of header\n";
    fwrite(cdata, sizeof(CHAR), (strlen(cdata) + 1), plotdata_ptr);

    /**************************************/
    /******* Write the data section *******/
    /**************************************/

#if defined(IPG_TGT)

    /*
     * Write the data section in the IPG Legacy format...
     */

    end = (LONG) rsp_free_wavemem_p;
    base = (LONG) rsp_start_wavemem_p;

    plotdata_hdr.wavegen_size = end - base;
    plotdata_hdr.wavegen_addr = rsp_start_wavemem_p;

    for (i = 0; i < WF_MAX_PROCESSORS; i++)
    {
        plotdata_hdr.instr_addr[i] =
            (long *) &rsp_start_instmem_p[4 * rsp_instmem_start[i]];
        plotdata_hdr.instr_size[i] = (LONG) (8 * (rsp_instmem_next[i] -
                                                  rsp_instmem_start[i]));
    }

    /* Write the data packet with the data packet set */
    bytes = fwrite((PDIpgPlotDataHeader *) &plotdata_hdr, 
                   (INT) 1,
                   (INT) sizeof(PDIpgPlotDataHeader),
                   (FILE *) plotdata_ptr);

    /* 
     * So that this routine will also work on the hardware; copy the data 
     * from hardware memory to local memory where it can be written as 
     * bytes.
     */
    wmem = (SHORT *) WTAlloc(end - base + 1);
    if (wmem == NULL)
    {
        WriteError("Insufficient memory for normal program operation.");
        WTFree(filename);
        WTFree(psdname);
        CloseFile(plotdata_ptr);
        return FAILURE;
    }
    for (jj = 0; jj < ((end - base) / 2); jj++)
        wmem[jj] = rsp_start_wavemem_p[jj];

    /* 
     * Write the waveform memory to the PD/PDX file.
     */
    bytes = fwrite((SHORT *) wmem, 
                   (INT) 1,
                   (INT) (end - base),
                   (FILE *) plotdata_ptr);
    WTFree(wmem);			/* Free the node */

    /* 
     * Traverse through the Hardware Queues, copy the information to 
     * memory, and write it to the PD/PDX file.
     */
    for (i = 0; i < WF_MAX_PROCESSORS; i++)
    {
        LONG nlen;
        SHORT *base;

        nlen = plotdata_hdr.instr_size[i] / sizeof(SHORT);
        imem = (SHORT *) WTAlloc(plotdata_hdr.instr_size[i] + 1);
        if (imem == NULL)
        {
            WriteError("Insufficient memory for normal program operation.");
            WTFree(filename);
            WTFree(psdname);
            CloseFile(plotdata_ptr);
            return FAILURE;
        }
        base = (SHORT *) &rsp_start_instmem_p[4 * rsp_instmem_start[i]];
        for (jj = 0; jj < nlen; jj++)
            imem[jj] = base[jj];
        bytes = fwrite((SHORT *) imem,
                       (INT) 1,
                       (INT) plotdata_hdr.instr_size[i],
                       (FILE *) plotdata_ptr);
        WTFree(imem);
    }

#elif defined(MGD_TGT)

    /*
     * Write the data section in the LxMGD 1.0 format...
     */

    fwrite(psdname, sizeof(CHAR), strlen(psdname), plotdata_ptr);
    fwrite("\n", sizeof(CHAR), 2, plotdata_ptr);	/* append \n\0 */
    cdata = "Created by PDStoreData()\n";
    fwrite(cdata, sizeof(CHAR), (strlen(cdata) + 1), plotdata_ptr);

    sprintf(cbuf, "n16-size = %d\n", sizeof(n16));
    fwrite(cbuf, sizeof(CHAR), (strlen(cbuf) + 1), plotdata_ptr);
    sprintf(cbuf, "n32-size = %d\n", sizeof(n32));
    fwrite(cbuf, sizeof(CHAR), (strlen(cbuf) + 1), plotdata_ptr);
    sprintf(cbuf, "pointer-size = %d\n", sizeof(void *));
    fwrite(cbuf, sizeof(CHAR), (strlen(cbuf) + 1), plotdata_ptr);

    n32_datum = NUM_WAVEGEN;
    fwrite(&n32_datum, sizeof(n32_datum), 1, plotdata_ptr);
    n32_datum = NUM_DSP;
    fwrite(&n32_datum, sizeof(n32_datum), 1, plotdata_ptr);

    for (i = 0; i < NUM_WAVEGEN; i++)
    {
        cdata = gen_name[i];
        fwrite(cdata, sizeof(CHAR), strlen(cdata), plotdata_ptr);
        fwrite("\n", sizeof(CHAR), 2, plotdata_ptr);    /* append \n\0 */
        n32_datum = i;
        fwrite(&n32_datum, sizeof(n32_datum), 1, plotdata_ptr);
        n32_datum = seq_dsp_index[i];
        fwrite(&n32_datum, sizeof(n32_datum), 1, plotdata_ptr);
    }

    for (i = 0; i < NUM_DSP; i++)
    {
        pointer_datum = (void *) (((char *) rsp_p->wavemem[i].base_ptr) -
                                  ((char *) rsp_p->shadow_mem_base_ptr[i]));
        fwrite(&pointer_datum, sizeof(pointer_datum), 1, plotdata_ptr);
        n32_datum = rsp_p->wavemem[i].wp_count;
        fwrite(&n32_datum, sizeof(n32_datum), 1, plotdata_ptr);
        n32_datum = rsp_p->wavemem[i].max_wp;
        fwrite(&n32_datum, sizeof(n32_datum), 1, plotdata_ptr);

        wmem = (n16 *) WTAlloc(sizeof(n16) * rsp_p->wavemem[i].wp_count);
        if (wmem == NULL)
        {
            WriteError("Insufficient memory for normal program operation.");
            WTFree(filename);
            WTFree(psdname);
            CloseFile(plotdata_ptr);
            return FAILURE;
        }
        for (jj = 0; jj < rsp_p->wavemem[i].wp_count; jj++)
            wmem[jj] = rsp_p->wavemem[i].base_ptr[jj];
        fwrite(wmem, sizeof(n16), rsp_p->wavemem[i].wp_count, plotdata_ptr);
        WTFree(wmem);
    }

    for (i = 0; i < NUM_WAVEGEN; i++)
    {
        pointer_datum = (void *) (rsp_p->instmem[i].base_ptr);
        fwrite(&pointer_datum, sizeof(pointer_datum), 1, plotdata_ptr);
        n32_datum = rsp_p->instmem[i].inst_count;
        fwrite(&n32_datum, sizeof(n32_datum), 1, plotdata_ptr);
        n32_datum = rsp_p->instmem[i].max_inst;
        fwrite(&n32_datum, sizeof(n32_datum), 1, plotdata_ptr);

        imem = (SeqInstruction0_s *) WTAlloc(sizeof(SeqInstruction0_s) * 
                                             rsp_p->instmem[i].inst_count);
        if (imem == NULL)
        {
            WriteError("Insufficient memory for normal program operation.");
            WTFree(filename);
            WTFree(psdname);
            CloseFile(plotdata_ptr);
            return FAILURE;
        }
        for (jj = 0; jj < rsp_p->instmem[i].inst_count; jj++)
        {
            imem[jj].period    = rsp_p->instmem[i].base_ptr[jj].period;
            imem[jj].waveform  = rsp_p->instmem[i].base_ptr[jj].waveform;
            imem[jj].amplitude = rsp_p->instmem[i].base_ptr[jj].amplitude;
        }
        fwrite(imem, sizeof(SeqInstruction0_s), 
               rsp_p->instmem[i].inst_count, plotdata_ptr);
        WTFree(imem);
    }

#endif

    /**************************************/
    /******* Close the PD/PDX file ********/
    /**************************************/

    WTFree(filename);
    WTFree(psdname);
    if (CloseFile(plotdata_ptr) != SUCCESS)
        return FAILURE;

    return SUCCESS;
}   /* end PDStoreData() */

/*****************************************************************************/
/* 
 * Construct the PD filename from the rootfile...
 */

static CHAR *
constructPDFileName( CHAR *rootfile )
{
    CHAR *newname;

    if (rootfile == NULL)
    {
        WriteError("PSD target program name missing.");
        return NULL;
    }

    newname = (CHAR *) WTAlloc(strlen(rootfile) + 4);
    if (newname == NULL)
    {
        WriteError("Insufficient memory for normal program operation.");
        return NULL;
    }
    strcpy(newname, rootfile);
    strcat(newname,".pd");

    return newname;
}   /* end constructPDFileName() */

/*****************************************************************************/
/* 
 * Construct the PSD Name from the rootfile...
 */

static CHAR *
constructPSDName( CHAR *rootfile )
{
    CHAR *newname, *tail;

    if (rootfile == NULL)
    {
        WriteError("PSD target program name missing.");
        return NULL;
    }

    tail = strstr(rootfile, ".ipg");
    if (tail == NULL)
        tail = strstr(rootfile, ".mgd");
    if (tail == NULL)
        tail = strstr(rootfile, ".tgt");

    if (tail == NULL)
    {
        WriteError("Incompatible PSD target program name.");
        return NULL;
    }
    newname = (CHAR *) WTAlloc((tail - rootfile) + 1);
    if (newname == NULL)
    {
        WriteError("Insufficient memory for normal program operation.");
        return NULL;
    }
    strncpy(newname, rootfile, (tail - rootfile));
    newname[(tail - rootfile)] = '\0';

    return newname;
}   /* end constructPSDName() */

