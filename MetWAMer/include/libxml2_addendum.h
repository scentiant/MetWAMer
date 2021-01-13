/* libxml2_addendum.h
 * Michael Sparks (mespar1@gmail.com)
 * Last modified : 20 July 2013
 *
 * Various macro definitions that might have an equivalent in
 * the libxml2 API already, but I was unable to find them.
 *
 * Copyright (c) 2007,2013 Michael E Sparks
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

#ifndef LIBXML2_ADDENDUM_H
#define LIBXML2_ADDENDUM_H

#include "errors.h"

/* Adjusts pointer argument to its child named "target", *
 * which is expected to exist.                           */
#define DESCEND_TO_CHILD(myxmlNodePtr,target) { \
  myxmlNodePtr=myxmlNodePtr->children; \
  while(xmlStrcmp(myxmlNodePtr->name,(const xmlChar*)target)) \
    if((myxmlNodePtr=myxmlNodePtr->next)==NULL) { \
      (void)snprintf(errbuff,MAXERRLEN, \
        "Err (): Can't descend to child \"%s\"!\n",target); \
      FATALERROR(errbuff) \
    } \
}

/* Adjusts pointer argument to its child named "target", *
 * which may or may not exist.                           */
#define ATTEMPT_DESCENT_TO_CHILD(myxmlNodePtr,target) { \
  myxmlNodePtr=myxmlNodePtr->children; \
  while(myxmlNodePtr!=NULL&& \
        xmlStrcmp(myxmlNodePtr->name,(const xmlChar*)target)) \
    myxmlNodePtr=myxmlNodePtr->next; \
}

/* Adjusts pointer argument to its parent */
#define ASCEND_TO_PARENT(myxmlNodePtr) { \
  myxmlNodePtr=myxmlNodePtr->parent; \
}

/* Finds sibling with name "target" */
#define FIND_SIBLING(myxmlNodePtr,target) { \
  myxmlNodePtr=myxmlNodePtr->next; \
  while(xmlStrcmp(myxmlNodePtr->name,(const xmlChar*)target)) \
    if((myxmlNodePtr=myxmlNodePtr->next)==NULL) \
      break; \
}

/* Makes expression of char-valued attribute recording clearer */
#define RECORD_ATTR_VAL_AS_CHAR(myxmlNodePtr,dest,attr) { \
  xmlChar *gen_attr=NULL; \
  if((gen_attr=xmlGetProp(myxmlNodePtr,(const xmlChar*)attr))==NULL) \
    FATALERROR("Err (): Attribute DNE!)\n") \
  dest=(char)gen_attr[0]; \
  xmlFree(gen_attr); \
}

/* Makes expression of char-string--valued attribute recording clearer */
#define RECORD_ATTR_VAL_AS_CHAR_STRING(myxmlNodePtr,dest,attr) { \
  xmlChar *gen_attr=NULL; \
  int attrvallen=-1; \
  if((gen_attr=xmlGetProp(myxmlNodePtr,(const xmlChar*)attr))==NULL) \
    FATALERROR("Err (): Attribute DNE!)\n") \
  else \
    attrvallen=strlen((const char*)gen_attr); \
  if((dest=(char*)malloc(sizeof(char)*(attrvallen+1)))==NULL) \
    FATALERROR("Err (): Insufficient Memory!)\n") \
  else \
    (void)strcpy(dest,(const char*)gen_attr); \
  xmlFree(gen_attr); \
}

/* Makes expression of int-valued attribute recording clearer */
#define RECORD_ATTR_VAL_AS_INT(myxmlNodePtr,dest,attr) { \
  xmlChar *gen_attr=NULL; \
  if((gen_attr=xmlGetProp(myxmlNodePtr,(const xmlChar*)attr))==NULL) \
    FATALERROR("Err (): Attribute DNE!)\n") \
  dest=atoi((const char*)gen_attr); \
  xmlFree(gen_attr); \
}

/* Add attribute attrname with int-valued attrval to node myxmlNodePtr */
#define SET_INT_AS_ATTR_VAL(myxmlNodePtr,attrname,attrval) { \
  (void)snprintf(errbuff,MAXERRLEN,"%i",attrval); \
  if(xmlNewProp(myxmlNodePtr,BAD_CAST attrname,BAD_CAST errbuff)==NULL) \
    FATALERROR("Err (): Unable to create a new attribute!\n") \
}
  
/* Add attribute attrname with real-valued attrval to node myxmlNodePtr */
#define SET_REAL_AS_ATTR_VAL(myxmlNodePtr,attrname,attrval) { \
  (void)snprintf(errbuff,MAXERRLEN,"%f",attrval); \
  if(xmlNewProp(myxmlNodePtr,BAD_CAST attrname,BAD_CAST errbuff)==NULL) \
    FATALERROR("Err (): Unable to create a new attribute!\n") \
}

/* Add attribute attrname with charstring- *
 * valued attrval to node myxmlNodePtr.    */
#define SET_STRING_AS_ATTR_VAL(myxmlNodePtr,attrname,attrval) { \
  if(xmlNewProp(myxmlNodePtr,BAD_CAST attrname,BAD_CAST attrval)==NULL) \
    FATALERROR("Err (): Unable to create a new attribute!\n") \
}
  
#endif
