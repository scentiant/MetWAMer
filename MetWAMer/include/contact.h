/* contact.h
 * Michael E Sparks (mespar1@gmail.com)
 * Last modified: 20 August 2007
 *
 * How to find me.
 *
 * Copyright (c) 2007 Michael E Sparks
 * All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#ifndef CONTACT_H
#define CONTACT_H

#include "errors.h"

char maintainer[]="Michael E Sparks",
     email[]="mespar1@gmail.com",
     url[]="http://brendelgroup.org/mespar1";

#define BUGGERS_OFF { \
  (void)snprintf(errbuff,MAXERRLEN,"\n  Please direct bug reports to:\n\n\
    %s (%s)\n\
    %s\n\n",maintainer,email,url); \
  NONFATALERROR(errbuff) \
}

#endif
