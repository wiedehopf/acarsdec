/*
 *  Copyright (c) 2007,2025 Thierry Leconte
 *
 *   
 *   This code is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU Library General Public License version 2
 *   published by the Free Software Foundation.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Library General Public License for more details.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "acarsdec.h"
#include "acars.h"

#define CEILING(x,y) (((x) + (y) - 1) / (y))

#define MSKFREQMARK 2400
#define MSKFREQSPACE 1200
#define MSKFREQCNTR  ((MSKFREQSPACE+MSKFREQMARK)/2)

#define BITLEN CEILING(INTRATE, MSKFREQSPACE)
#define MFLTOVER 240U
#define MFLTLEN (BITLEN * MFLTOVER + 1)

#if INTRATE % MSKFREQSPACE
 #warning INTRATE is not a multiple of MSQFREQSPACE, code may give odd results
#endif

static float h[MFLTLEN];

int initMsk(channel_t *ch)
{
	unsigned int i;

	ch->MskClk = ch->MskS = 0;
	ch->MskPhi = ch->MskDphi = ch->MskDf = 0;

	ch->idx = 0;
	ch->inb = calloc(BITLEN, sizeof(*ch->inb));
	if (ch->inb == NULL) {
		perror(NULL);
		return -1;
	}

	if (ch->chn == 0) {
		/* precompute half-wave matched filter table */
		for (i = 0; i < MFLTLEN; i++)
			h[i] = sin(M_PI * MSKFREQSPACE * (float)i  / INTRATE / MFLTOVER) ;
	}

	return 0;
}

// ACARS is LSb first
static inline void putbit(float v, channel_t *ch)
{
	ch->outbits >>= 1;
	if (v > 0)
		ch->outbits |= 0x80;

	if (--ch->nbits == 0)
		decodeAcars(ch);
}

static const float PLLKi = 71e-7/BITLEN;
static const float PLLKp = 60e-3/BITLEN;

void demodMSK(channel_t *ch, int len)
{
	/* MSK demod */
	int n;
	unsigned int idx = ch->idx;
	float p = ch->MskPhi;
	float s;

	for (n = 0; n < len; n++) {
		float in;
		float complex v;
		unsigned int j, o;

		s = 2.0 * M_PI * (float)(MSKFREQCNTR) / INTRATE + ch->MskDphi;

		/* bit clock */
		ch->MskClk += s;
		if (ch->MskClk > 3 * M_PI / 2.0) {
			float dphi;
			float vo, lvl;

			ch->MskClk -= 3 * M_PI / 2.0;

			/* matched filter */
			o = MFLTOVER * (ch->MskClk / s );
			if (o > MFLTOVER)
				o = MFLTOVER;
			for (v = 0, j = 0; j < BITLEN; j++, o += MFLTOVER)
				v += h[o] * ch->inb[(j + idx) % BITLEN];

			/* normalize */
			lvl = cabsf(v) + 1e-8F;
			v /= lvl;

			/* update magnitude exp moving average. Average over last 8 bits */
			ch->MskMag = ch->MskMag - (1.0F/8.0F * (ch->MskMag - lvl));

			if (ch->MskS & 1) {
				// Q
				vo = cimagf(v);
				dphi = (vo >= 0) ? -crealf(v) : crealf(v);
			} else {
				// I
				vo = crealf(v);
				dphi = (vo >= 0) ? cimagf(v) : -cimagf(v);
			}

			putbit((ch->MskS & 2) ? -vo : vo, ch);
			ch->MskS++;

			// lock on signal once ACARS header has been / is being heard
			if (PREKEY != ch->Acarsstate || ch->count) {
				/* PLL as a PI controller */
				ch->MskDf += PLLKi * dphi;
				ch->MskDphi = ch->MskDf + PLLKp * dphi;
			}
			else	// otherwise don't even try to lock. XXX REVISIT: use a 2nd order / 1st order PLL split?
				ch->MskDf = ch->MskDphi = 0
		}

		/* VCO */
		p += s;
		if (p >= 2.0 * M_PI)
			p -= 2.0 * M_PI;

		/* mixer */
		in = ch->dm_buffer[n];
		ch->inb[idx] = in * cexp(-p * I);
		idx = (idx + 1) % BITLEN;
	}

	ch->idx = idx;
	ch->MskPhi = p;
}
