/**
 ****************************************************************************
 * <P> XML.c - implementation file for basic XML parser written in ANSI C++
 * for portability. It works by using recursion and a node tree for breaking
 * down the elements of an XML document.  </P>
 *
 * @version     V2.24
 * @author      Frank Vanden Berghen
 *
 * NOTE:
 *
 *   If you add "#define STRICT_PARSING", on the first line of this file
 *   the parser will see the following XML-stream:
 *      <a><b>some text</b><b>other text    </a>
 *   as an error. Otherwise, this tring will be equivalent to:
 *      <a><b>some text</b><b>other text</b></a>
 *
 * NOTE:
 *
 *   If you add "#define APPROXIMATE_PARSING" on the first line of this file
 *   the parser will see the following XML-stream:
 *     <data name="n1">
 *     <data name="n2">
 *     <data name="n3" />
 *   as equivalent to the following XML-stream:
 *     <data name="n1" />
 *     <data name="n2" />
 *     <data name="n3" />
 *   This can be useful for badly-formed XML-streams but prevent the use
 *   of the following XML-stream (problem is: tags at contiguous levels
 *   have the same names):
 *     <data name="n1">
 *        <data name="n2">
 *            <data name="n3" />
 *        </data>
 *     </data>
 *
 * NOTE:
 *
 *   If you add "#define _XMLPARSER_NO_MESSAGEBOX_" on the first line of this file
 *   the "openFileHelper" function will always display error messages inside the
 *   console instead of inside a message-box-window. Message-box-windows are
 *   available on windows 9x/NT/2000/XP/Vista only.
 *
 * BSD license:
 * Copyright (c) 2002, Frank Vanden Berghen
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Frank Vanden Berghen nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************
 */
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include "xmlParser.h"
#ifdef _XMLWINDOWS
//#ifdef _DEBUG
//#define _CRTDBG_MAP_ALLOC
//#include <crtdbg.h>
//#endif
#define WIN32_LEAN_AND_MEAN
#include <Windows.h> // to have IsTextUnicode, MultiByteToWideChar, WideCharToMultiByte to handle unicode files
                     // to have "MessageBoxA" to display error messages for openFilHelper
#endif

#include <memory.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <wchar.h>

XMLCSTR XMLNode::getVersion() { return _T("v2.23"); }
void free_XMLDLL(void *t){free(t);}

static char strictUTF8Parsing=1, guessUnicodeChars=1, dropWhiteSpace=1;

inline int mmin( const int t1, const int t2 ) { return t1 < t2 ? t1 : t2; }

// You can modify the initialization of the variable "XMLClearTags" below
// to change the clearTags that are currently recognized by the library.
// The number on the second columns is the length of the string inside the
// first column. The "<!DOCTYPE" declaration must be the second in the list.
static ALLXMLClearTag XMLClearTags[] =
{
    {    _T("<![CDATA["),9,  _T("]]>")      },
    {    _T("<!DOCTYPE"),9,  _T(">")        },
    {    _T("<PRE>")    ,5,  _T("</PRE>")   },
    {    _T("<Script>") ,8,  _T("</Script>")},
    {    _T("<!--")     ,4,  _T("-->")      },
    {    NULL           ,0,  NULL           }
};
ALLXMLClearTag* XMLNode::getClearTagTable() { return XMLClearTags; }

// You can modify the initialization of the variable "XMLEntities" below
// to change the character entities that are currently recognized by the library.
// The number on the second columns is the length of the string inside the
// first column. Additionally, the syntaxes "&#xA0;" and "&#160;" are recognized.
typedef struct { XMLCSTR s; int l; XMLCHAR c;} XMLCharacterEntity;
static XMLCharacterEntity XMLEntities[] =
{
    { _T("&amp;" ), 5, _T('&' )},
    { _T("&lt;"  ), 4, _T('<' )},
    { _T("&gt;"  ), 4, _T('>' )},
    { _T("&quot;"), 6, _T('\"')},
    { _T("&apos;"), 6, _T('\'')},
    { NULL        , 0, '\0'    }
};

// When rendering the XMLNode to a string (using the "createXMLString" function),
// you can ask for a beautiful formatting. This formatting is using the
// following indentation character:
#define INDENTCHAR _T('\t')

// The following function parses the XML errors into a user friendly string.
// You can edit this to change the output language of the library to something else.
XMLCSTR XMLNode::getError(XMLError xerror)
{
    switch (xerror)
    {
    case eXMLErrorNone:                  return _T("No error");
    case eXMLErrorMissingEndTag:         return _T("Warning: Unmatched end tag");
    case eXMLErrorEmpty:                 return _T("Error: No XML data");
    case eXMLErrorFirstNotStartTag:      return _T("Error: First token not start tag");
    case eXMLErrorMissingTagName:        return _T("Error: Missing start tag name");
    case eXMLErrorMissingEndTagName:     return _T("Error: Missing end tag name");
    case eXMLErrorNoMatchingQuote:       return _T("Error: Unmatched quote");
    case eXMLErrorUnmatchedEndTag:       return _T("Error: Unmatched end tag");
    case eXMLErrorUnmatchedEndClearTag:  return _T("Error: Unmatched clear tag end");
    case eXMLErrorUnexpectedToken:       return _T("Error: Unexpected token found");
    case eXMLErrorInvalidTag:            return _T("Error: Invalid tag found");
    case eXMLErrorNoElements:            return _T("Error: No elements found");
    case eXMLErrorFileNotFound:          return _T("Error: File not found");
    case eXMLErrorFirstTagNotFound:      return _T("Error: First Tag not found");
    case eXMLErrorUnknownCharacterEntity:return _T("Error: Unknown character entity");
    case eXMLErrorCharConversionError:   return _T("Error: unable to convert between UNICODE and MultiByte chars");
    case eXMLErrorCannotOpenWriteFile:   return _T("Error: unable to open file for writing");
    case eXMLErrorCannotWriteFile:       return _T("Error: cannot write into file");

    case eXMLErrorBase64DataSizeIsNotMultipleOf4: return _T("Warning: Base64-string length is not a multiple of 4");
    case eXMLErrorBase64DecodeTruncatedData:      return _T("Warning: Base64-string is truncated");
    case eXMLErrorBase64DecodeIllegalCharacter:   return _T("Error: Base64-string contains an illegal character");
    case eXMLErrorBase64DecodeBufferTooSmall:     return _T("Error: Base64 decode output buffer is too small");
    };
    return _T("Unknown");
}

// Here is an abstraction layer to access some common string manipulation functions.
// The abstraction layer is currently working for gcc, Microsoft Visual Studio 6.0,
// Microsoft Visual Studio .NET, CC (sun compiler) and Borland C++.
// If you plan to "port" the library to a new system/compiler, all you have to do is
// to edit the following lines.
#ifdef XML_NO_WIDE_CHAR
char myIsTextUnicode(const void *b, int len) { return FALSE; }
#else
    #if defined (UNDER_CE) || !defined(WIN32)
    char myIsTextUnicode(const void *b, int len) // inspired by the Wine API: RtlIsTextUnicode
    {
#ifdef sun
        // for SPARC processors: wchar_t* buffers must always be alligned, otherwise it's a char* buffer.
        if ((((unsigned long)b)%sizeof(wchar_t))!=0) return FALSE;
#endif
        const wchar_t *s=(const wchar_t*)b;

        // buffer too small:
        if (len<(int)sizeof(wchar_t)) return FALSE;

        // odd length test
        if (len&1) return FALSE;

        /* only checks the first 256 characters */
        len=mmin(256,len/sizeof(wchar_t));

        // Check for the special byte order:
        if (*s == 0xFFFE) return FALSE;     // IS_TEXT_UNICODE_REVERSE_SIGNATURE;
        if (*s == 0xFEFF) return TRUE;      // IS_TEXT_UNICODE_SIGNATURE

        // checks for ASCII characters in the UNICODE stream
        int i,stats=0;
        for (i=0; i<len; i++) if (s[i]<=(unsigned short)255) stats++;
        if (stats>len/2) return TRUE;

        // Check for UNICODE NULL chars
        for (i=0; i<len; i++) if (!s[i]) return TRUE;

        return FALSE;
    }
    #else
    char myIsTextUnicode(const void *b,int l) { return (char)IsTextUnicode((CONST LPVOID)b,l,NULL); };
    #endif
#endif

#ifdef _XMLWINDOWS
// for Microsoft Visual Studio 6.0 and Microsoft Visual Studio .NET,
    #ifdef _XMLUNICODE
        wchar_t *myMultiByteToWideChar(const char *s,int l)
        {
            int i;
            if (strictUTF8Parsing)  i=(int)MultiByteToWideChar(CP_UTF8,0             ,s,l,NULL,0);
            else                    i=(int)MultiByteToWideChar(CP_ACP ,MB_PRECOMPOSED,s,l,NULL,0);
            if (i<0) return NULL;
            wchar_t *d=(wchar_t *)malloc((i+1)*sizeof(XMLCHAR));
            if (strictUTF8Parsing)  i=(int)MultiByteToWideChar(CP_UTF8,0             ,s,l,d,i);
            else                    i=(int)MultiByteToWideChar(CP_ACP ,MB_PRECOMPOSED,s,l,d,i);
            d[i]=0;
            return d;
        }
    #else
        char *myWideCharToMultiByte(const wchar_t *s,int l)
        {
            UINT codePage=CP_ACP; if (strictUTF8Parsing) codePage=CP_UTF8;
            int i=(int)WideCharToMultiByte(codePage,  // code page
                0,                       // performance and mapping flags
                s,                       // wide-character string
                l,                       // number of chars in string
                NULL,                       // buffer for new string
                0,                       // size of buffer
                NULL,                    // default for unmappable chars
                NULL                     // set when default char used
                );
            if (i<0) return NULL;
            char *d=(char*)malloc(i+1);
            WideCharToMultiByte(codePage,  // code page
                0,                       // performance and mapping flags
                s,                       // wide-character string
                l,                       // number of chars in string
                d,                       // buffer for new string
                i,                       // size of buffer
                NULL,                    // default for unmappable chars
                NULL                     // set when default char used
                );
            d[i]=0;
            return d;
        }
    #endif
    #ifdef __BORLANDC__
    int _strnicmp(char *c1, char *c2, int l){ return strnicmp(c1,c2,l);}
    #endif
#else
// for gcc and CC
    #ifdef XML_NO_WIDE_CHAR
        char *myWideCharToMultiByte(const wchar_t *s, int l) { return NULL; }
    #else
        char *myWideCharToMultiByte(const wchar_t *s, int l)
        {
            const wchar_t *ss=s;
            int i=(int)wcsrtombs(NULL,&ss,0,NULL);
            if (i<0) return NULL;
            char *d=(char *)malloc(i+1);
            wcsrtombs(d,&s,i,NULL);
            d[i]=0;
            return d;
        }
    #endif
    #ifdef _XMLUNICODE
        wchar_t *myMultiByteToWideChar(const char *s, int l)
        {
            const char *ss=s;
            int i=(int)mbsrtowcs(NULL,&ss,0,NULL);
            if (i<0) return NULL;
            wchar_t *d=(wchar_t *)malloc((i+1)*sizeof(wchar_t));
            mbsrtowcs(d,&s,l,NULL);
            d[i]=0;
            return d;
        }
        int _tcslen(XMLCSTR c)   { return wcslen(c); }
        #ifdef sun
        // for CC
           #include <widec.h>
           int _tcsnicmp(XMLCSTR c1, XMLCSTR c2, int l) { return wsncasecmp(c1,c2,l);}
           int _tcsicmp(XMLCSTR c1, XMLCSTR c2) { return wscasecmp(c1,c2); }
        #else
        // for gcc
           int _tcsnicmp(XMLCSTR c1, XMLCSTR c2, int l) { return wcsncasecmp(c1,c2,l);}
           int _tcsicmp(XMLCSTR c1, XMLCSTR c2) { return wcscasecmp(c1,c2); }
        #endif
        XMLSTR _tcsstr(XMLCSTR c1, XMLCSTR c2) { return (XMLSTR)wcsstr(c1,c2); }
        XMLSTR _tcscpy(XMLSTR c1, XMLCSTR c2) { return (XMLSTR)wcscpy(c1,c2); }
        FILE *_tfopen(XMLCSTR filename,XMLCSTR mode)
        {
            char *filenameAscii=myWideCharToMultiByte(filename,0);
            FILE *f;
            if (mode[0]==_T('r')) f=fopen(filenameAscii,"rb");
            else                  f=fopen(filenameAscii,"wb");
            free(filenameAscii);
            return f;
        }
    #else
        FILE *_tfopen(XMLCSTR filename,XMLCSTR mode) { return fopen(filename,mode); }
        int _tcslen(XMLCSTR c)   { return strlen(c); }
        int _tcsnicmp(XMLCSTR c1, XMLCSTR c2, int l) { return strncasecmp(c1,c2,l);}
        int _tcsicmp(XMLCSTR c1, XMLCSTR c2) { return strcasecmp(c1,c2); }
        XMLSTR _tcsstr(XMLCSTR c1, XMLCSTR c2) { return (XMLSTR)strstr(c1,c2); }
        XMLSTR _tcscpy(XMLSTR c1, XMLCSTR c2) { return (XMLSTR)strcpy(c1,c2); }
    #endif
    int _strnicmp(const char *c1,const char *c2, int l) { return strncasecmp(c1,c2,l);}
#endif

/////////////////////////////////////////////////////////////////////////
//      Here start the core implementation of the XMLParser library    //
/////////////////////////////////////////////////////////////////////////

// You should normally not change anything below this point.
// For your own information, I suggest that you read the openFileHelper below:
XMLNode XMLNode::openFileHelper(XMLCSTR filename, XMLCSTR tag)
{
    // guess the value of the global parameter "strictUTF8Parsing"
    // (the guess is based on the first 200 bytes of the file).
    FILE *f=_tfopen(filename,_T("rb"));
    if (f)
    {
        char bb[205];
        int l=(int)fread(bb,1,200,f);
        setGlobalOptions(guessUnicodeChars,guessUTF8ParsingParameterValue(bb,l),dropWhiteSpace);
        fclose(f);
    }

    // parse the file
    XMLResults pResults;
    XMLNode xnode=XMLNode::parseFile(filename,tag,&pResults);

    // display error message (if any)
    if (pResults.error != eXMLErrorNone)
    {
        // create message
        char message[2000],*s1=(char*)"",*s3=(char*)""; XMLCSTR s2=_T("");
        if (pResults.error==eXMLErrorFirstTagNotFound) { s1=(char*)"First Tag should be '"; s2=tag; s3=(char*)"'.\n"; }
        sprintf(message,
#ifdef _XMLUNICODE
            "XML Parsing error inside file '%S'.\n%S\nAt line %i, column %i.\n%s%S%s"
#else
            "XML Parsing error inside file '%s'.\n%s\nAt line %i, column %i.\n%s%s%s"
#endif
            ,filename,XMLNode::getError(pResults.error),pResults.nLine,pResults.nColumn,s1,s2,s3);

        // display message
#if defined(WIN32) && !defined(UNDER_CE) && !defined(_XMLPARSER_NO_MESSAGEBOX_)
        MessageBoxA(NULL,message,"XML Parsing error",MB_OK|MB_ICONERROR|MB_TOPMOST);
#else
        printf("%s",message);
#endif
        exit(255);
    }
    return xnode;
}

#ifndef _XMLUNICODE
// If "strictUTF8Parsing=0" then we assume that all characters have the same length of 1 byte.
// If "strictUTF8Parsing=1" then the characters have different lengths (from 1 byte to 4 bytes).
// This table is used as lookup-table to know the length of a character (in byte) based on the
// content of the first byte of the character.
// (note: if you modify this, you must always have XML_utf8ByteTable[0]=0 ).
static const char XML_utf8ByteTable[256] =
{
    //  0 1 2 3 4 5 6 7 8 9 a b c d e f
    0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x00
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x10
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x20
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x30
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x40
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x50
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x60
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x70End of ASCII range
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x80 0x80 to 0xc1 invalid
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0x90
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0xa0
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,// 0xb0
    1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,// 0xc0 0xc2 to 0xdf 2 byte
    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,// 0xd0
    3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,// 0xe0 0xe0 to 0xef 3 byte
    4,4,4,4,4,1,1,1,1,1,1,1,1,1,1,1 // 0xf0 0xf0 to 0xf4 4 byte, 0xf5 and higher invalid
};
static const char XML_asciiByteTable[256] =
{
    0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
};
static const char *XML_ByteTable=(const char *)XML_utf8ByteTable; // the default is "strictUTF8Parsing=1"
#endif

XMLError XMLNode::writeToFile(XMLCSTR filename, const char *encoding, char nFormat) const
{
    if (!d) return eXMLErrorNone;
    int i;
    XMLSTR t=createXMLString(nFormat,&i);
    FILE *f=_tfopen(filename,_T("wb"));
    if (!f) return eXMLErrorCannotOpenWriteFile;
#ifdef _XMLUNICODE
    unsigned char h[2]={ 0xFF, 0xFE };
    if (!fwrite(h,2,1,f)) return eXMLErrorCannotWriteFile;
    if ((!isDeclaration())&&((d->lpszName)||(!getChildNode().isDeclaration())))
    {
        if (!fwrite(_T("<?xml version=\"1.0\" encoding=\"utf-16\"?>\n"),sizeof(wchar_t)*40,1,f))
            return eXMLErrorCannotWriteFile;
    }
#else
    if ((!isDeclaration())&&((d->lpszName)||(!getChildNode().isDeclaration())))
    {
        if ((!encoding)||(XML_ByteTable==XML_utf8ByteTable))
        {
            // header so that windows recognize the file as UTF-8:
            unsigned char h[3]={0xEF,0xBB,0xBF};
            if (!fwrite(h,3,1,f)) return eXMLErrorCannotWriteFile;
            if (!fwrite("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n",39,1,f)) return eXMLErrorCannotWriteFile;
        }
        else
            if (fprintf(f,"<?xml version=\"1.0\" encoding=\"%s\"?>\n",encoding)<0) return eXMLErrorCannotWriteFile;
    } else
    {
        if (XML_ByteTable==XML_utf8ByteTable) // test if strictUTF8Parsing==1"
        {
            unsigned char h[3]={0xEF,0xBB,0xBF}; if (!fwrite(h,3,1,f)) return eXMLErrorCannotWriteFile;
        }
    }
#endif
    if (!fwrite(t,sizeof(XMLCHAR)*i,1,f)) return eXMLErrorCannotWriteFile;
    if (fclose(f)!=0) return eXMLErrorCannotWriteFile;
    free(t);
    return eXMLErrorNone;
}

// Duplicate a given string.
XMLSTR stringDup(XMLCSTR lpszData, int cbData)
{
    if (lpszData==NULL) return NULL;

    XMLSTR lpszNew;
    if (cbData==0) cbData=(int)_tcslen(lpszData);
    lpszNew = (XMLSTR)malloc((cbData+1) * sizeof(XMLCHAR));
    if (lpszNew)
    {
        memcpy(lpszNew, lpszData, (cbData) * sizeof(XMLCHAR));
        lpszNew[cbData] = (XMLCHAR)NULL;
    }
    return lpszNew;
}

XMLNode XMLNode::emptyXMLNode;
XMLClear XMLNode::emptyXMLClear={ NULL, NULL, NULL};
XMLAttribute XMLNode::emptyXMLAttribute={ NULL, NULL};

// Enumeration used to decipher what type a token is
typedef enum XMLTokenTypeTag
{
    eTokenText = 0,
    eTokenQuotedText,
    eTokenTagStart,         /* "<"            */
    eTokenTagEnd,           /* "</"           */
    eTokenCloseTag,         /* ">"            */
    eTokenEquals,           /* "="            */
    eTokenDeclaration,      /* "<?"           */
    eTokenShortHandClose,   /* "/>"           */
    eTokenClear,
    eTokenError
} XMLTokenType;

// Main structure used for parsing XML
typedef struct XML
{
    XMLCSTR                lpXML;
    XMLCSTR                lpszText;
    int                    nIndex,nIndexMissigEndTag;
    enum XMLError          error;
    XMLCSTR                lpEndTag;
    int                    cbEndTag;
    XMLCSTR                lpNewElement;
    int                    cbNewElement;
    int                    nFirst;
} XML;

typedef struct
{
    ALLXMLClearTag *pClr;
    XMLCSTR     pStr;
} NextToken;

// Enumeration used when parsing attributes
typedef enum Attrib
{
    eAttribName = 0,
    eAttribEquals,
    eAttribValue
} Attrib;

// Enumeration used when parsing elements to dictate whether we are currently
// inside a tag
typedef enum Status
{
    eInsideTag = 0,
    eOutsideTag
} Status;

// private (used while rendering):
XMLSTR toXMLString(XMLSTR dest,XMLCSTR source)
{
    XMLSTR dd=dest;
    XMLCHAR ch;
    XMLCharacterEntity *entity;
    while ((ch=*source))
    {
        entity=XMLEntities;
        do
        {
            if (ch==entity->c) {_tcscpy(dest,entity->s); dest+=entity->l; source++; goto out_of_loop1; }
            entity++;
        } while(entity->s);
#ifdef _XMLUNICODE
        *(dest++)=*(source++);
#else
        switch(XML_ByteTable[(unsigned char)ch])
        {
        case 4: *(dest++)=*(source++);
        case 3: *(dest++)=*(source++);
        case 2: *(dest++)=*(source++);
        case 1: *(dest++)=*(source++);
        }
#endif
out_of_loop1:
        ;
    }
    *dest=0;
    return dd;
}

// private (used while rendering):
int lengthXMLString(XMLCSTR source)
{
    int r=0;
    XMLCharacterEntity *entity;
    XMLCHAR ch;
    while ((ch=*source))
    {
        entity=XMLEntities;
        do
        {
            if (ch==entity->c) { r+=entity->l; source++; goto out_of_loop1; }
            entity++;
        } while(entity->s);
#ifdef _XMLUNICODE
        r++; source++;
#else
        ch=XML_ByteTable[(unsigned char)ch]; r+=ch; source+=ch;
#endif
out_of_loop1:
        ;
    }
    return r;
}

XMLSTR toXMLString(XMLCSTR source)
{
    XMLSTR dest=(XMLSTR)malloc((lengthXMLString(source)+1)*sizeof(XMLCHAR));
    return toXMLString(dest,source);
}

XMLSTR toXMLStringFast(XMLSTR *dest,int *destSz, XMLCSTR source)
{
    int l=lengthXMLString(source)+1;
    if (l>*destSz) { *destSz=l; *dest=(XMLSTR)realloc(*dest,l*sizeof(XMLCHAR)); }
    return toXMLString(*dest,source);
}

// private:
XMLSTR fromXMLString(XMLCSTR s, int lo, XML *pXML)
{
    // This function is the opposite of the function "toXMLString". It decodes the escape
    // sequences &amp;, &quot;, &apos;, &lt;, &gt; and replace them by the characters
    // &,",',<,>. This function is used internally by the XML Parser. All the calls to
    // the XML library will always gives you back "decoded" strings.
    //
    // in: string (s) and length (lo) of string
    // out:  new allocated string converted from xml
    if (!s) return NULL;

    int ll=0,j;
    XMLSTR d;
    XMLCSTR ss=s;
    XMLCharacterEntity *entity;
    while ((lo>0)&&(*s))
    {
        if (*s==_T('&'))
        {
            if ((lo>2)&&(s[1]==_T('#')))
            {
                s+=2; lo-=2;
                if ((*s==_T('X'))||(*s==_T('x'))) { s++; lo--; }
                while ((*s)&&(*s!=_T(';'))&&((lo--)>0)) s++;
                if (*s!=_T(';'))
                {
                    pXML->error=eXMLErrorUnknownCharacterEntity;
                    return NULL;
                }
                s++; lo--;
            } else
            {
                entity=XMLEntities;
                do
                {
                    if ((lo>=entity->l)&&(_tcsnicmp(s,entity->s,entity->l)==0)) { s+=entity->l; lo-=entity->l; break; }
                    entity++;
                } while(entity->s);
                if (!entity->s)
                {
                    pXML->error=eXMLErrorUnknownCharacterEntity;
                    return NULL;
                }
            }
        } else
        {
#ifdef _XMLUNICODE
            s++; lo--;
#else
            j=XML_ByteTable[(unsigned char)*s]; s+=j; lo-=j; ll+=j-1;
#endif
        }
        ll++;
    }

    d=(XMLSTR)malloc((ll+1)*sizeof(XMLCHAR));
    s=d;
    while (ll-->0)
    {
        if (*ss==_T('&'))
        {
            if (ss[1]==_T('#'))
            {
                ss+=2; j=0;
                if ((*ss==_T('X'))||(*ss==_T('x')))
                {
                    ss++;
                    while (*ss!=_T(';'))
                    {
                        if ((*ss>=_T('0'))&&(*ss<=_T('9'))) j=(j<<4)+*ss-_T('0');
                        else if ((*ss>=_T('A'))&&(*ss<=_T('F'))) j=(j<<4)+*ss-_T('A')+10;
                        else if ((*ss>=_T('a'))&&(*ss<=_T('f'))) j=(j<<4)+*ss-_T('a')+10;
                        else { free((void*)s); pXML->error=eXMLErrorUnknownCharacterEntity;return NULL;}
                        ss++;
                    }
                } else
                {
                    while (*ss!=_T(';'))
                    {
                        if ((*ss>=_T('0'))&&(*ss<=_T('9'))) j=(j*10)+*ss-_T('0');
                        else { free((void*)s); pXML->error=eXMLErrorUnknownCharacterEntity;return NULL;}
                        ss++;
                    }
                }
                (*d++)=(XMLCHAR)j; ss++;
            } else
            {
                entity=XMLEntities;
                do
                {
                    if (_tcsnicmp(ss,entity->s,entity->l)==0) { *(d++)=entity->c; ss+=entity->l; break; }
                    entity++;
                } while(entity->s);
            }
        } else
        {
#ifdef _XMLUNICODE
            *(d++)=*(ss++);
#else
            switch(XML_ByteTable[(unsigned char)*ss])
            {
            case 4: *(d++)=*(ss++); ll--;
            case 3: *(d++)=*(ss++); ll--;
            case 2: *(d++)=*(ss++); ll--;
            case 1: *(d++)=*(ss++);
            }
#endif
        }
    }
    *d=0;
    return (XMLSTR)s;
}

#define XML_isSPACECHAR(ch) ((ch==_T('\n'))||(ch==_T(' '))||(ch== _T('\t'))||(ch==_T('\r')))

// private:
char myTagCompare(XMLCSTR cclose, XMLCSTR copen)
// !!!! WARNING strange convention&:
// return 0 if equals
// return 1 if different
{
    if (!cclose) return 1;
    int l=(int)_tcslen(cclose);
    if (_tcsnicmp(cclose, copen, l)!=0) return 1;
    const XMLCHAR c=copen[l];
    if (XML_isSPACECHAR(c)||
        (c==_T('/' ))||
        (c==_T('<' ))||
        (c==_T('>' ))||
        (c==_T('=' ))) return 0;
    return 1;
}

// Obtain the next character from the string.
static inline XMLCHAR getNextChar(XML *pXML)
{
    XMLCHAR ch = pXML->lpXML[pXML->nIndex];
#ifdef _XMLUNICODE
    if (ch!=0) pXML->nIndex++;
#else
    pXML->nIndex+=XML_ByteTable[(unsigned char)ch];
#endif
    return ch;
}

// Find the next token in a string.
// pcbToken contains the number of characters that have been read.
static NextToken GetNextToken(XML *pXML, int *pcbToken, enum XMLTokenTypeTag *pType)
{
    NextToken        result;
    XMLCHAR            ch;
    XMLCHAR            chTemp;
    int              indexStart,nFoundMatch,nIsText=FALSE;
    result.pClr=NULL; // prevent warning

    // Find next non-white space character
    do { indexStart=pXML->nIndex; ch=getNextChar(pXML); } while XML_isSPACECHAR(ch);

    if (ch)
    {
        // Cache the current string pointer
        result.pStr = &pXML->lpXML[indexStart];

        // First check whether the token is in the clear tag list (meaning it
        // does not need formatting).
        ALLXMLClearTag *ctag=XMLClearTags;
        do
        {
            if (_tcsnicmp(ctag->lpszOpen, result.pStr, ctag->openTagLen)==0)
            {
                result.pClr=ctag;
                pXML->nIndex+=ctag->openTagLen-1;
                *pType=eTokenClear;
                return result;
            }
            ctag++;
        } while(ctag->lpszOpen);

        // If we didn't find a clear tag then check for standard tokens
        switch(ch)
        {
        // Check for quotes
        case _T('\''):
        case _T('\"'):
            // Type of token
            *pType = eTokenQuotedText;
            chTemp = ch;

            // Set the size
            nFoundMatch = FALSE;

            // Search through the string to find a matching quote
            while((ch = getNextChar(pXML)))
            {
                if (ch==chTemp) { nFoundMatch = TRUE; break; }
                if (ch==_T('<')) break;
            }

            // If we failed to find a matching quote
            if (nFoundMatch == FALSE)
            {
                pXML->nIndex=indexStart+1;
                nIsText=TRUE;
                break;
            }

//  4.02.2002
//            if (FindNonWhiteSpace(pXML)) pXML->nIndex--;

            break;

        // Equals (used with attribute values)
        case _T('='):
            *pType = eTokenEquals;
            break;

        // Close tag
        case _T('>'):
            *pType = eTokenCloseTag;
            break;

        // Check for tag start and tag end
        case _T('<'):

            // Peek at the next character to see if we have an end tag '</',
            // or an xml declaration '<?'
            chTemp = pXML->lpXML[pXML->nIndex];

            // If we have a tag end...
            if (chTemp == _T('/'))
            {
                // Set the type and ensure we point at the next character
                getNextChar(pXML);
                *pType = eTokenTagEnd;
            }

            // If we have an XML declaration tag
            else if (chTemp == _T('?'))
            {

                // Set the type and ensure we point at the next character
                getNextChar(pXML);
                *pType = eTokenDeclaration;
            }

            // Otherwise we must have a start tag
            else
            {
                *pType = eTokenTagStart;
            }
            break;

        // Check to see if we have a short hand type end tag ('/>').
        case _T('/'):

            // Peek at the next character to see if we have a short end tag '/>'
            chTemp = pXML->lpXML[pXML->nIndex];

            // If we have a short hand end tag...
            if (chTemp == _T('>'))
            {
                // Set the type and ensure we point at the next character
                getNextChar(pXML);
                *pType = eTokenShortHandClose;
                break;
            }

            // If we haven't found a short hand closing tag then drop into the
            // text process

        // Other characters
        default:
            nIsText = TRUE;
        }

        // If this is a TEXT node
        if (nIsText)
        {
            // Indicate we are dealing with text
            *pType = eTokenText;
            while((ch = getNextChar(pXML)))
            {
                if XML_isSPACECHAR(ch)
                {
                    indexStart++; break;

                } else if (ch==_T('/'))
                {
                    // If we find a slash then this maybe text or a short hand end tag
                    // Peek at the next character to see it we have short hand end tag
                    ch=pXML->lpXML[pXML->nIndex];
                    // If we found a short hand end tag then we need to exit the loop
                    if (ch==_T('>')) { pXML->nIndex--; break; }

                } else if ((ch==_T('<'))||(ch==_T('>'))||(ch==_T('=')))
                {
                    pXML->nIndex--; break;
                }
            }
        }
        *pcbToken = pXML->nIndex-indexStart;
    } else
    {
        // If we failed to obtain a valid character
        *pcbToken = 0;
        *pType = eTokenError;
        result.pStr=NULL;
    }

    return result;
}

XMLCSTR XMLNode::updateName_WOSD(XMLCSTR lpszName)
{
    if (d->lpszName&&(lpszName!=d->lpszName)) free((void*)d->lpszName);
    d->lpszName=lpszName;
    return lpszName;
}

// private:
XMLNode::XMLNode(struct XMLNodeDataTag *p){ d=p; (p->ref_count)++; }
XMLNode::XMLNode(XMLNodeData *pParent, XMLCSTR lpszName, char isDeclaration)
{
    d=(XMLNodeData*)malloc(sizeof(XMLNodeData));
    d->ref_count=1;

    d->lpszName=NULL;
    d->nChild= 0;
    d->nText = 0;
    d->nClear = 0;
    d->nAttribute = 0;

    d->isDeclaration = isDeclaration;

    d->pParent = pParent;
    d->pChild= NULL;
    d->pText= NULL;
    d->pClear= NULL;
    d->pAttribute= NULL;
    d->pOrder= NULL;

    updateName_WOSD(lpszName);
}

XMLNode XMLNode::createXMLTopNode_WOSD(XMLCSTR lpszName, char isDeclaration) { return XMLNode(NULL,lpszName,isDeclaration); }
XMLNode XMLNode::createXMLTopNode(XMLCSTR lpszName, char isDeclaration) { return XMLNode(NULL,stringDup(lpszName),isDeclaration); }

#define MEMORYINCREASE 50

static inline void *myRealloc(void *p, int newsize, int memInc, int sizeofElem)
{
    if (p==NULL) { if (memInc) return malloc(memInc*sizeofElem); return malloc(sizeofElem); }
    if ((memInc==0)||((newsize%memInc)==0)) p=realloc(p,(newsize+memInc)*sizeofElem);
//    if (!p)
//    {
//        printf("XMLParser Error: Not enough memory! Aborting...\n"); exit(220);
//    }
    return p;
}

// private:
int XMLNode::findPosition(XMLNodeData *d, int index, XMLElementType xtype)
{
    if (index<0) return -1;
    int i=0,j=(int)((index<<2)+xtype),*o=d->pOrder; while (o[i]!=j) i++; return i;
}

// private:
// update "order" information when deleting a content of a XMLNode
int XMLNode::removeOrderElement(XMLNodeData *d, XMLElementType t, int index)
{
    int n=d->nChild+d->nText+d->nClear, *o=d->pOrder,i=findPosition(d,index,t);
    memmove(o+i, o+i+1, (n-i)*sizeof(int));
    for (;i<n;i++)
        if ((o[i]&3)==(int)t) o[i]-=4;
    // We should normally do:
    // d->pOrder=(int)realloc(d->pOrder,n*sizeof(int));
    // but we skip reallocation because it's too time consuming.
    // Anyway, at the end, it will be free'd completely at once.
    return i;
}

void *XMLNode::addToOrder(int memoryIncrease,int *_pos, int nc, void *p, int size, XMLElementType xtype)
{
    //  in: *_pos is the position inside d->pOrder ("-1" means "EndOf")
    // out: *_pos is the index inside p
    p=myRealloc(p,(nc+1),memoryIncrease,size);
    int n=d->nChild+d->nText+d->nClear;
    d->pOrder=(int*)myRealloc(d->pOrder,n+1,memoryIncrease*3,sizeof(int));
    int pos=*_pos,*o=d->pOrder;

    if ((pos<0)||(pos>=n)) { *_pos=nc; o[n]=(int)((nc<<2)+xtype); return p; }

    int i=pos;
    memmove(o+i+1, o+i, (n-i)*sizeof(int));

    while ((pos<n)&&((o[pos]&3)!=(int)xtype)) pos++;
    if (pos==n) { *_pos=nc; o[n]=(int)((nc<<2)+xtype); return p; }

    o[i]=o[pos];
    for (i=pos+1;i<=n;i++) if ((o[i]&3)==(int)xtype) o[i]+=4;

    *_pos=pos=o[pos]>>2;
    memmove(((char*)p)+(pos+1)*size,((char*)p)+pos*size,(nc-pos)*size);

    return p;
}

// Add a child node to the given element.
XMLNode XMLNode::addChild_priv(int memoryIncrease, XMLCSTR lpszName, char isDeclaration, int pos)
{
    if (!lpszName) return emptyXMLNode;
    d->pChild=(XMLNode*)addToOrder(memoryIncrease,&pos,d->nChild,d->pChild,sizeof(XMLNode),eNodeChild);
    d->pChild[pos].d=NULL;
    d->pChild[pos]=XMLNode(d,lpszName,isDeclaration);
    d->nChild++;
    return d->pChild[pos];
}

// Add an attribute to an element.
XMLAttribute *XMLNode::addAttribute_priv(int memoryIncrease,XMLCSTR lpszName, XMLCSTR lpszValuev)
{
    if (!lpszName) return &emptyXMLAttribute;
    int nc=d->nAttribute;
    d->pAttribute=(XMLAttribute*)myRealloc(d->pAttribute,(nc+1),memoryIncrease,sizeof(XMLAttribute));
    XMLAttribute *pAttr=d->pAttribute+nc;
    pAttr->lpszName = lpszName;
    pAttr->lpszValue = lpszValuev;
    d->nAttribute++;
    return pAttr;
}

// Add text to the element.
XMLCSTR XMLNode::addText_priv(int memoryIncrease, XMLCSTR lpszValue, int pos)
{
    if (!lpszValue) return NULL;
    d->pText=(XMLCSTR*)addToOrder(memoryIncrease,&pos,d->nText,d->pText,sizeof(XMLSTR),eNodeText);
    d->pText[pos]=lpszValue;
    d->nText++;
    return lpszValue;
}

// Add clear (unformatted) text to the element.
XMLClear *XMLNode::addClear_priv(int memoryIncrease, XMLCSTR lpszValue, XMLCSTR lpszOpen, XMLCSTR lpszClose, int pos)
{
    if (!lpszValue) return &emptyXMLClear;
    d->pClear=(XMLClear *)addToOrder(memoryIncrease,&pos,d->nClear,d->pClear,sizeof(XMLClear),eNodeClear);
    XMLClear *pNewClear=d->pClear+pos;
    pNewClear->lpszValue = lpszValue;
    if (!lpszOpen) lpszOpen=getClearTagTable()->lpszOpen;
    if (!lpszClose) lpszOpen=getClearTagTable()->lpszClose;
    pNewClear->lpszOpenTag = lpszOpen;
    pNewClear->lpszCloseTag = lpszClose;
    d->nClear++;
    return pNewClear;
}

// private:
// Parse a clear (unformatted) type node.
char XMLNode::parseClearTag(void *px, ALLXMLClearTag *pClear)
{
    XML *pXML=(XML *)px;
    int cbTemp=0;
    XMLCSTR lpszTemp=NULL;
    XMLCSTR lpXML=&pXML->lpXML[pXML->nIndex];
    static XMLCSTR docTypeEnd=_T("]>");

    // Find the closing tag
    // Seems the <!DOCTYPE need a better treatment so lets handle it
    if (pClear->lpszOpen==XMLClearTags[1].lpszOpen)
    {
        XMLCSTR pCh=lpXML;
        while (*pCh)
        {
            if (*pCh==_T('<')) { pClear->lpszClose=docTypeEnd; lpszTemp=_tcsstr(lpXML,docTypeEnd); break; }
            else if (*pCh==_T('>')) { lpszTemp=pCh; break; }
#ifdef _XMLUNICODE
            pCh++;
#else
            pCh+=XML_ByteTable[(unsigned char)(*pCh)];
#endif
        }
    } else lpszTemp=_tcsstr(lpXML, pClear->lpszClose);

    if (lpszTemp)
    {
        // Cache the size and increment the index
        cbTemp = (int)(lpszTemp - lpXML);

        pXML->nIndex += cbTemp+(int)_tcslen(pClear->lpszClose);

        // Add the clear node to the current element
        addClear_priv(MEMORYINCREASE,stringDup(lpXML,cbTemp), pClear->lpszOpen, pClear->lpszClose,-1);
        return 0;
    }

    // If we failed to find the end tag
    pXML->error = eXMLErrorUnmatchedEndClearTag;
    return 1;
}

void XMLNode::exactMemory(XMLNodeData *d)
{
    if (d->pOrder)     d->pOrder=(int*)realloc(d->pOrder,(d->nChild+d->nText+d->nClear)*sizeof(int));
    if (d->pChild)     d->pChild=(XMLNode*)realloc(d->pChild,d->nChild*sizeof(XMLNode));
    if (d->pAttribute) d->pAttribute=(XMLAttribute*)realloc(d->pAttribute,d->nAttribute*sizeof(XMLAttribute));
    if (d->pText)      d->pText=(XMLCSTR*)realloc(d->pText,d->nText*sizeof(XMLSTR));
    if (d->pClear)     d->pClear=(XMLClear *)realloc(d->pClear,d->nClear*sizeof(XMLClear));
}

char XMLNode::maybeAddTxT(void *pa, XMLCSTR tokenPStr)
{
    XML *pXML=(XML *)pa;
    XMLCSTR lpszText=pXML->lpszText;
    if (!lpszText) return 0;
    if (dropWhiteSpace) while (XML_isSPACECHAR(*lpszText)&&(lpszText!=tokenPStr)) lpszText++;
    int cbText = (int)(tokenPStr - lpszText);
    if (!cbText) { pXML->lpszText=NULL; return 0; }
    if (dropWhiteSpace) { cbText--; while ((cbText)&&XML_isSPACECHAR(lpszText[cbText])) cbText--; cbText++; }
    if (!cbText) { pXML->lpszText=NULL; return 0; }
    lpszText=fromXMLString(lpszText,cbText,pXML);
    if (!lpszText) return 1;
    addText_priv(MEMORYINCREASE,lpszText,-1);
    pXML->lpszText=NULL;
    return 0;
}
// private:
// Recursively parse an XML element.
int XMLNode::ParseXMLElement(void *pa)
{
    XML *pXML=(XML *)pa;
    int cbToken;
    enum XMLTokenTypeTag type;
    NextToken token;
    XMLCSTR lpszTemp=NULL;
    int cbTemp=0;
    char nDeclaration;
    XMLNode pNew;
    enum Status status; // inside or outside a tag
    enum Attrib attrib = eAttribName;

    assert(pXML);

    // If this is the first call to the function
    if (pXML->nFirst)
    {
        // Assume we are outside of a tag definition
        pXML->nFirst = FALSE;
        status = eOutsideTag;
    } else
    {
        // If this is not the first call then we should only be called when inside a tag.
        status = eInsideTag;
    }

    // Iterate through the tokens in the document
    for(;;)
    {
        // Obtain the next token
        token = GetNextToken(pXML, &cbToken, &type);

        if (type != eTokenError)
        {
            // Check the current status
            switch(status)
            {

            // If we are outside of a tag definition
            case eOutsideTag:

                // Check what type of token we obtained
                switch(type)
                {
                // If we have found text or quoted text
                case eTokenText:
                case eTokenCloseTag:          /* '>'         */
                case eTokenShortHandClose:    /* '/>'        */
                case eTokenQuotedText:
                case eTokenEquals:
                    break;

                // If we found a start tag '<' and declarations '<?'
                case eTokenTagStart:
                case eTokenDeclaration:

                    // Cache whether this new element is a declaration or not
                    nDeclaration = (type == eTokenDeclaration);

                    // If we have node text then add this to the element
                    if (maybeAddTxT(pXML,token.pStr)) return FALSE;

                    // Find the name of the tag
                    token = GetNextToken(pXML, &cbToken, &type);

                    // Return an error if we couldn't obtain the next token or
                    // it wasnt text
                    if (type != eTokenText)
                    {
                        pXML->error = eXMLErrorMissingTagName;
                        return FALSE;
                    }

                    // If we found a new element which is the same as this
                    // element then we need to pass this back to the caller..

#ifdef APPROXIMATE_PARSING
                    if (d->lpszName &&
                        myTagCompare(d->lpszName, token.pStr) == 0)
                    {
                        // Indicate to the caller that it needs to create a
                        // new element.
                        pXML->lpNewElement = token.pStr;
                        pXML->cbNewElement = cbToken;
                        return TRUE;
                    } else
#endif
                    {
                        // If the name of the new element differs from the name of
                        // the current element we need to add the new element to
                        // the current one and recurse
                        pNew = addChild_priv(MEMORYINCREASE,stringDup(token.pStr,cbToken), nDeclaration,-1);

                        while (!pNew.isEmpty())
                        {
                            // Callself to process the new node.  If we return
                            // FALSE this means we dont have any more
                            // processing to do...

                            if (!pNew.ParseXMLElement(pXML)) return FALSE;
                            else
                            {
                                // If the call to recurse this function
                                // evented in a end tag specified in XML then
                                // we need to unwind the calls to this
                                // function until we find the appropriate node
                                // (the element name and end tag name must
                                // match)
                                if (pXML->cbEndTag)
                                {
                                    // If we are back at the root node then we
                                    // have an unmatched end tag
                                    if (!d->lpszName)
                                    {
                                        pXML->error=eXMLErrorUnmatchedEndTag;
                                        return FALSE;
                                    }

                                    // If the end tag matches the name of this
                                    // element then we only need to unwind
                                    // once more...

                                    if (myTagCompare(d->lpszName, pXML->lpEndTag)==0)
                                    {
                                        pXML->cbEndTag = 0;
                                    }

                                    return TRUE;
                                } else
                                    if (pXML->cbNewElement)
                                    {
                                        // If the call indicated a new element is to
                                        // be created on THIS element.

                                        // If the name of this element matches the
                                        // name of the element we need to create
                                        // then we need to return to the caller
                                        // and let it process the element.

                                        if (myTagCompare(d->lpszName, pXML->lpNewElement)==0)
                                        {
                                            return TRUE;
                                        }

                                        // Add the new element and recurse
                                        pNew = addChild_priv(MEMORYINCREASE,stringDup(pXML->lpNewElement,pXML->cbNewElement),0,-1);
                                        pXML->cbNewElement = 0;
                                    }
                                    else
                                    {
                                        // If we didn't have a new element to create
                                        pNew = emptyXMLNode;

                                    }
                            }
                        }
                    }
                    break;

                // If we found an end tag
                case eTokenTagEnd:

                    // If we have node text then add this to the element
                    if (maybeAddTxT(pXML,token.pStr)) return FALSE;

                    // Find the name of the end tag
                    token = GetNextToken(pXML, &cbTemp, &type);

                    // The end tag should be text
                    if (type != eTokenText)
                    {
                        pXML->error = eXMLErrorMissingEndTagName;
                        return FALSE;
                    }
                    lpszTemp = token.pStr;

                    // After the end tag we should find a closing tag
                    token = GetNextToken(pXML, &cbToken, &type);
                    if (type != eTokenCloseTag)
                    {
                        pXML->error = eXMLErrorMissingEndTagName;
                        return FALSE;
                    }
                    pXML->lpszText=pXML->lpXML+pXML->nIndex;

                    // We need to return to the previous caller.  If the name
                    // of the tag cannot be found we need to keep returning to
                    // caller until we find a match
                    if (myTagCompare(d->lpszName, lpszTemp) != 0)
#ifdef STRICT_PARSING
                    {
                        pXML->error=eXMLErrorUnmatchedEndTag;
                        pXML->nIndexMissigEndTag=pXML->nIndex;
                        return FALSE;
                    }
#else
                    {
                        pXML->error=eXMLErrorMissingEndTag;
                        pXML->nIndexMissigEndTag=pXML->nIndex;
                        pXML->lpEndTag = lpszTemp;
                        pXML->cbEndTag = cbTemp;
                    }
#endif

                    // Return to the caller
                    exactMemory(d);
                    return TRUE;

                // If we found a clear (unformatted) token
                case eTokenClear:
                    // If we have node text then add this to the element
                    if (maybeAddTxT(pXML,token.pStr)) return FALSE;
                    if (parseClearTag(pXML, token.pClr)) return FALSE;
                    pXML->lpszText=pXML->lpXML+pXML->nIndex;
                    break;

                default:
                    break;
                }
                break;

            // If we are inside a tag definition we need to search for attributes
            case eInsideTag:

                // Check what part of the attribute (name, equals, value) we
                // are looking for.
                switch(attrib)
                {
                // If we are looking for a new attribute
                case eAttribName:

                    // Check what the current token type is
                    switch(type)
                    {
                    // If the current type is text...
                    // Eg.  'attribute'
                    case eTokenText:
                        // Cache the token then indicate that we are next to
                        // look for the equals
                        lpszTemp = token.pStr;
                        cbTemp = cbToken;
                        attrib = eAttribEquals;
                        break;

                    // If we found a closing tag...
                    // Eg.  '>'
                    case eTokenCloseTag:
                        // We are now outside the tag
                        status = eOutsideTag;
                        pXML->lpszText=pXML->lpXML+pXML->nIndex;
                        break;

                    // If we found a short hand '/>' closing tag then we can
                    // return to the caller
                    case eTokenShortHandClose:
                        exactMemory(d);
                        pXML->lpszText=pXML->lpXML+pXML->nIndex;
                        return TRUE;

                    // Errors...
                    case eTokenQuotedText:    /* '"SomeText"'   */
                    case eTokenTagStart:      /* '<'            */
                    case eTokenTagEnd:        /* '</'           */
                    case eTokenEquals:        /* '='            */
                    case eTokenDeclaration:   /* '<?'           */
                    case eTokenClear:
                        pXML->error = eXMLErrorUnexpectedToken;
                        return FALSE;
                    default: break;
                    }
                    break;

                // If we are looking for an equals
                case eAttribEquals:
                    // Check what the current token type is
                    switch(type)
                    {
                    // If the current type is text...
                    // Eg.  'Attribute AnotherAttribute'
                    case eTokenText:
                        // Add the unvalued attribute to the list
                        addAttribute_priv(MEMORYINCREASE,stringDup(lpszTemp,cbTemp), NULL);
                        // Cache the token then indicate.  We are next to
                        // look for the equals attribute
                        lpszTemp = token.pStr;
                        cbTemp = cbToken;
                        break;

                    // If we found a closing tag 'Attribute >' or a short hand
                    // closing tag 'Attribute />'
                    case eTokenShortHandClose:
                    case eTokenCloseTag:
                        // If we are a declaration element '<?' then we need
                        // to remove extra closing '?' if it exists
                        pXML->lpszText=pXML->lpXML+pXML->nIndex;

                        if (d->isDeclaration &&
                            (lpszTemp[cbTemp-1]) == _T('?'))
                        {
                            cbTemp--;
                        }

                        if (cbTemp)
                        {
                            // Add the unvalued attribute to the list
                            addAttribute_priv(MEMORYINCREASE,stringDup(lpszTemp,cbTemp), NULL);
                        }

                        // If this is the end of the tag then return to the caller
                        if (type == eTokenShortHandClose)
                        {
                            exactMemory(d);
                            return TRUE;
                        }

                        // We are now outside the tag
                        status = eOutsideTag;
                        break;

                    // If we found the equals token...
                    // Eg.  'Attribute ='
                    case eTokenEquals:
                        // Indicate that we next need to search for the value
                        // for the attribute
                        attrib = eAttribValue;
                        break;

                    // Errors...
                    case eTokenQuotedText:    /* 'Attribute "InvalidAttr"'*/
                    case eTokenTagStart:      /* 'Attribute <'            */
                    case eTokenTagEnd:        /* 'Attribute </'           */
                    case eTokenDeclaration:   /* 'Attribute <?'           */
                    case eTokenClear:
                        pXML->error = eXMLErrorUnexpectedToken;
                        return FALSE;
                    default: break;
                    }
                    break;

                // If we are looking for an attribute value
                case eAttribValue:
                    // Check what the current token type is
                    switch(type)
                    {
                    // If the current type is text or quoted text...
                    // Eg.  'Attribute = "Value"' or 'Attribute = Value' or
                    // 'Attribute = 'Value''.
                    case eTokenText:
                    case eTokenQuotedText:
                        // If we are a declaration element '<?' then we need
                        // to remove extra closing '?' if it exists
                        if (d->isDeclaration &&
                            (token.pStr[cbToken-1]) == _T('?'))
                        {
                            cbToken--;
                        }

                        if (cbTemp)
                        {
                            // Add the valued attribute to the list
                            if (type==eTokenQuotedText) { token.pStr++; cbToken-=2; }
                            XMLCSTR attrVal=token.pStr;
                            if (attrVal)
                            {
                                attrVal=fromXMLString(attrVal,cbToken,pXML);
                                if (!attrVal) return FALSE;
                            }
                            addAttribute_priv(MEMORYINCREASE,stringDup(lpszTemp,cbTemp),attrVal);
                        }

                        // Indicate we are searching for a new attribute
                        attrib = eAttribName;
                        break;

                    // Errors...
                    case eTokenTagStart:        /* 'Attr = <'          */
                    case eTokenTagEnd:          /* 'Attr = </'         */
                    case eTokenCloseTag:        /* 'Attr = >'          */
                    case eTokenShortHandClose:  /* "Attr = />"         */
                    case eTokenEquals:          /* 'Attr = ='          */
                    case eTokenDeclaration:     /* 'Attr = <?'         */
                    case eTokenClear:
                        pXML->error = eXMLErrorUnexpectedToken;
                        return FALSE;
                        break;
                    default: break;
                    }
                }
            }
        }
        // If we failed to obtain the next token
        else
        {
            if ((!d->isDeclaration)&&(d->pParent))
            {
#ifdef STRICT_PARSING
                pXML->error=eXMLErrorUnmatchedEndTag;
#else
                pXML->error=eXMLErrorMissingEndTag;
#endif
                pXML->nIndexMissigEndTag=pXML->nIndex;
            }
            return FALSE;
        }
    }
}

// Count the number of lines and columns in an XML string.
static void CountLinesAndColumns(XMLCSTR lpXML, int nUpto, XMLResults *pResults)
{
    XMLCHAR ch;
    assert(lpXML);
    assert(pResults);

    struct XML xml={ lpXML,lpXML, 0, 0, eXMLErrorNone, NULL, 0, NULL, 0, TRUE };

    pResults->nLine = 1;
    pResults->nColumn = 1;
    while (xml.nIndex<nUpto)
    {
        ch = getNextChar(&xml);
        if (ch != _T('\n')) pResults->nColumn++;
        else
        {
            pResults->nLine++;
            pResults->nColumn=1;
        }
    }
}

// Parse XML and return the root element.
XMLNode XMLNode::parseString(XMLCSTR lpszXML, XMLCSTR tag, XMLResults *pResults)
{
    if (!lpszXML)
    {
        if (pResults)
        {
            pResults->error=eXMLErrorNoElements;
            pResults->nLine=0;
            pResults->nColumn=0;
        }
        return emptyXMLNode;
    }

    XMLNode xnode(NULL,NULL,FALSE);
    struct XML xml={ lpszXML, lpszXML, 0, 0, eXMLErrorNone, NULL, 0, NULL, 0, TRUE };

    // Create header element
    xnode.ParseXMLElement(&xml);
    enum XMLError error = xml.error;
    if ((xnode.nChildNode()==1)&&(xnode.nElement()==1)) xnode=xnode.getChildNode(); // skip the empty node

    // If no error occurred
    if ((error==eXMLErrorNone)||(error==eXMLErrorMissingEndTag))
    {
        XMLCSTR name=xnode.getName();
        if (tag&&_tcslen(tag)&&((!name)||(_tcsicmp(xnode.getName(),tag))))
        {
            XMLNode nodeTmp;
            int i=0;
            while (i<xnode.nChildNode())
            {
                nodeTmp=xnode.getChildNode(i);
                if (_tcsicmp(nodeTmp.getName(),tag)==0) break;
                if (nodeTmp.isDeclaration()) { xnode=nodeTmp; i=0; } else i++;
            }
            if (i>=xnode.nChildNode())
            {
                if (pResults)
                {
                    pResults->error=eXMLErrorFirstTagNotFound;
                    pResults->nLine=0;
                    pResults->nColumn=0;
                }
                return emptyXMLNode;
            }
            xnode=nodeTmp;
        }
    } else
    {
        // Cleanup: this will destroy all the nodes
        xnode = emptyXMLNode;
    }


    // If we have been given somewhere to place results
    if (pResults)
    {
        pResults->error = error;

        // If we have an error
        if (error!=eXMLErrorNone)
        {
            if (error==eXMLErrorMissingEndTag) xml.nIndex=xml.nIndexMissigEndTag;
            // Find which line and column it starts on.
            CountLinesAndColumns(xml.lpXML, xml.nIndex, pResults);
        }
    }
    return xnode;
}

XMLNode XMLNode::parseFile(XMLCSTR filename, XMLCSTR tag, XMLResults *pResults)
{
    if (pResults) { pResults->nLine=0; pResults->nColumn=0; }
    FILE *f=_tfopen(filename,_T("rb"));
    if (f==NULL) { if (pResults) pResults->error=eXMLErrorFileNotFound; return emptyXMLNode; }
    fseek(f,0,SEEK_END);
    int l=ftell(f),headerSz=0;
    if (!l) { if (pResults) pResults->error=eXMLErrorEmpty; return emptyXMLNode; }
    fseek(f,0,SEEK_SET);
    unsigned char *buf=(unsigned char*)malloc(l+1);
    fread(buf,l,1,f);
    fclose(f);
    buf[l]=0;
#ifdef _XMLUNICODE
    if (guessUnicodeChars)
    {
        if (!myIsTextUnicode(buf,l))
        {
            if ((buf[0]==0xef)&&(buf[1]==0xbb)&&(buf[2]==0xbf)) headerSz=3;
            XMLSTR b2=myMultiByteToWideChar((const char*)(buf+headerSz),l-headerSz);
            free(buf); buf=(unsigned char*)b2; headerSz=0;
        } else
        {
            if ((buf[0]==0xef)&&(buf[1]==0xff)) headerSz=2;
            if ((buf[0]==0xff)&&(buf[1]==0xfe)) headerSz=2;
        }
    }
#else
    if (guessUnicodeChars)
    {
        if (myIsTextUnicode(buf,l))
        {
            l/=sizeof(wchar_t);
            if ((buf[0]==0xef)&&(buf[1]==0xff)) headerSz=2;
            if ((buf[0]==0xff)&&(buf[1]==0xfe)) headerSz=2;
            char *b2=myWideCharToMultiByte((const wchar_t*)(buf+headerSz),l-headerSz);
            free(buf); buf=(unsigned char*)b2; headerSz=0;
        } else
        {
            if ((buf[0]==0xef)&&(buf[1]==0xbb)&&(buf[2]==0xbf)) headerSz=3;
        }
    }
#endif

    if (!buf) { if (pResults) pResults->error=eXMLErrorCharConversionError; return emptyXMLNode; }
    XMLNode x=parseString((XMLSTR)(buf+headerSz),tag,pResults);
    free(buf);
    return x;
}

static inline void charmemset(XMLSTR dest,XMLCHAR c,int l) { while (l--) *(dest++)=c; }
// private:
// Creates an user friendly XML string from a given element with
// appropriate white space and carriage returns.
//
// This recurses through all subnodes then adds contents of the nodes to the
// string.
int XMLNode::CreateXMLStringR(XMLNodeData *pEntry, XMLSTR lpszMarker, int nFormat)
{
    int nResult = 0;
    int cb;
    int cbElement;
    int nChildFormat=-1;
    int nElementI=pEntry->nChild+pEntry->nText+pEntry->nClear;
    int i,j;

    assert(pEntry);

#define LENSTR(lpsz) (lpsz ? _tcslen(lpsz) : 0)

    // If the element has no name then assume this is the head node.
    cbElement = (int)LENSTR(pEntry->lpszName);

    if (cbElement)
    {
        // "<elementname "
        cb = nFormat == -1 ? 0 : nFormat;

        if (lpszMarker)
        {
            if (cb) charmemset(lpszMarker, INDENTCHAR, sizeof(XMLCHAR)*cb);
            nResult = cb;
            lpszMarker[nResult++]=_T('<');
            if (pEntry->isDeclaration) lpszMarker[nResult++]=_T('?');
            _tcscpy(&lpszMarker[nResult], pEntry->lpszName);
            nResult+=cbElement;
            lpszMarker[nResult++]=_T(' ');

        } else
        {
            nResult+=cbElement+2+cb;
            if (pEntry->isDeclaration) nResult++;
        }

        // Enumerate attributes and add them to the string
        XMLAttribute *pAttr=pEntry->pAttribute;
        for (i=0; i<pEntry->nAttribute; i++)
        {
            // "Attrib
            cb = (int)LENSTR(pAttr->lpszName);
            if (cb)
            {
                if (lpszMarker) _tcscpy(&lpszMarker[nResult], pAttr->lpszName);
                nResult += cb;
                // "Attrib=Value "
                if (pAttr->lpszValue)
                {
                    cb=(int)lengthXMLString(pAttr->lpszValue);
                    if (lpszMarker)
                    {
                        lpszMarker[nResult]=_T('=');
                        lpszMarker[nResult+1]=_T('"');
                        if (cb) toXMLString(&lpszMarker[nResult+2],pAttr->lpszValue);
                        lpszMarker[nResult+cb+2]=_T('"');
                    }
                    nResult+=cb+3;
                }
                if (lpszMarker) lpszMarker[nResult] = _T(' ');
                nResult++;
            }
            pAttr++;
        }

        if (pEntry->isDeclaration)
        {
            if (lpszMarker)
            {
                lpszMarker[nResult-1]=_T('?');
                lpszMarker[nResult]=_T('>');
            }
            nResult++;
            if (nFormat!=-1)
            {
                if (lpszMarker) lpszMarker[nResult]=_T('\n');
                nResult++;
            }
        } else
            // If there are child nodes we need to terminate the start tag
            if (nElementI)
            {
                if (lpszMarker) lpszMarker[nResult-1]=_T('>');
                if (nFormat!=-1)
                {
                    if (lpszMarker) lpszMarker[nResult]=_T('\n');
                    nResult++;
                }
            } else nResult--;
    }

    // Calculate the child format for when we recurse.  This is used to
    // determine the number of spaces used for prefixes.
    if (nFormat!=-1)
    {
        if (cbElement&&(!pEntry->isDeclaration)) nChildFormat=nFormat+1;
        else nChildFormat=nFormat;
    }

    // Enumerate through remaining children
    for (i=0; i<nElementI; i++)
    {
        j=pEntry->pOrder[i];
        switch((XMLElementType)(j&3))
        {
        // Text nodes
        case eNodeText:
            {
                // "Text"
                XMLCSTR pChild=pEntry->pText[j>>2];
                cb = (int)lengthXMLString(pChild);
                if (cb)
                {
                    if (nFormat!=-1)
                    {
                        if (lpszMarker)
                        {
                            charmemset(&lpszMarker[nResult],INDENTCHAR,sizeof(XMLCHAR)*(nFormat + 1));
                            toXMLString(&lpszMarker[nResult+nFormat+1],pChild);
                            lpszMarker[nResult+nFormat+1+cb]=_T('\n');
                        }
                        nResult+=cb+nFormat+2;
                    } else
                    {
                        if (lpszMarker) toXMLString(&lpszMarker[nResult], pChild);
                        nResult += cb;
                    }
                }
                break;
            }

        // Clear type nodes
        case eNodeClear:
            {
                XMLClear *pChild=pEntry->pClear+(j>>2);
                // "OpenTag"
                cb = (int)LENSTR(pChild->lpszOpenTag);
                if (cb)
                {
                    if (nFormat!=-1)
                    {
                        if (lpszMarker)
                        {
                            charmemset(&lpszMarker[nResult], INDENTCHAR, sizeof(XMLCHAR)*(nFormat + 1));
                            _tcscpy(&lpszMarker[nResult+nFormat+1], pChild->lpszOpenTag);
                        }
                        nResult+=cb+nFormat+1;
                    }
                    else
                    {
                        if (lpszMarker)_tcscpy(&lpszMarker[nResult], pChild->lpszOpenTag);
                        nResult += cb;
                    }
                }

                // "OpenTag Value"
                cb = (int)LENSTR(pChild->lpszValue);
                if (cb)
                {
                    if (lpszMarker) _tcscpy(&lpszMarker[nResult], pChild->lpszValue);
                    nResult += cb;
                }

                // "OpenTag Value CloseTag"
                cb = (int)LENSTR(pChild->lpszCloseTag);
                if (cb)
                {
                    if (lpszMarker) _tcscpy(&lpszMarker[nResult], pChild->lpszCloseTag);
                    nResult += cb;
                }

                if (nFormat!=-1)
                {
                    if (lpszMarker) lpszMarker[nResult] = _T('\n');
                    nResult++;
                }
                break;
            }

        // Element nodes
        case eNodeChild:
            {
                // Recursively add child nodes
                nResult += CreateXMLStringR(pEntry->pChild[j>>2].d, lpszMarker ? lpszMarker + nResult : 0, nChildFormat);
                break;
            }
        default: break;
        }
    }

    if ((cbElement)&&(!pEntry->isDeclaration))
    {
        // If we have child entries we need to use long XML notation for
        // closing the element - "<elementname>blah blah blah</elementname>"
        if (nElementI)
        {
            // "</elementname>\0"
            if (lpszMarker)
            {
                if (nFormat != -1)
                {
                    if (nFormat)
                    {
                        charmemset(&lpszMarker[nResult], INDENTCHAR,sizeof(XMLCHAR)*nFormat);
                        nResult+=nFormat;
                    }
                }

                _tcscpy(&lpszMarker[nResult], _T("</"));
                nResult += 2;
                _tcscpy(&lpszMarker[nResult], pEntry->lpszName);
                nResult += cbElement;

                if (nFormat == -1)
                {
                    _tcscpy(&lpszMarker[nResult], _T(">"));
                    nResult++;
                } else
                {
                    _tcscpy(&lpszMarker[nResult], _T(">\n"));
                    nResult+=2;
                }
            } else
            {
                if (nFormat != -1) nResult+=cbElement+4+nFormat;
                else nResult+=cbElement+3;
            }
        } else
        {
            // If there are no children we can use shorthand XML notation -
            // "<elementname/>"
            // "/>\0"
            if (lpszMarker)
            {
                if (nFormat == -1)
                {
                    _tcscpy(&lpszMarker[nResult], _T("/>"));
                    nResult += 2;
                }
                else
                {
                    _tcscpy(&lpszMarker[nResult], _T("/>\n"));
                    nResult += 3;
                }
            }
            else
            {
                nResult += nFormat == -1 ? 2 : 3;
            }
        }
    }

    return nResult;
}

#undef LENSTR

// Create an XML string
// @param       int nFormat             - 0 if no formatting is required
//                                        otherwise nonzero for formatted text
//                                        with carriage returns and indentation.
// @param       int *pnSize             - [out] pointer to the size of the
//                                        returned string not including the
//                                        NULL terminator.
// @return      XMLSTR                  - Allocated XML string, you must free
//                                        this with free().
XMLSTR XMLNode::createXMLString(int nFormat, int *pnSize) const
{
    if (!d) { if (pnSize) *pnSize=0; return NULL; }

    XMLSTR lpszResult = NULL;
    int cbStr;

    // Recursively Calculate the size of the XML string
    if (!dropWhiteSpace) nFormat=0;
    nFormat = nFormat ? 0 : -1;
    cbStr = CreateXMLStringR(d, 0, nFormat);
    assert(cbStr);
    // Alllocate memory for the XML string + the NULL terminator and
    // create the recursively XML string.
    lpszResult=(XMLSTR)malloc((cbStr+1)*sizeof(XMLCHAR));
    CreateXMLStringR(d, lpszResult, nFormat);
    if (pnSize) *pnSize = cbStr;
    return lpszResult;
}

XMLNode::~XMLNode() { deleteNodeContent(); }

int XMLNode::detachFromParent(XMLNodeData *d)
{
    XMLNode *pa=d->pParent->pChild;
    int i=0;
    while (((void*)(pa[i].d))!=((void*)d)) i++;
    d->pParent->nChild--;
    if (d->pParent->nChild) memmove(pa+i,pa+i+1,(d->pParent->nChild-i)*sizeof(XMLNode));
    else { free(pa); d->pParent->pChild=NULL; }
    return removeOrderElement(d->pParent,eNodeChild,i);
}

void XMLNode::deleteNodeContent(char force)
{
    if (!d) return;
    (d->ref_count) --;
    if ((d->ref_count==0)||force)
    {
        int i;
        if (d->pParent) detachFromParent(d);
        for(i=0; i<d->nChild; i++) { d->pChild[i].d->pParent=NULL; d->pChild[i].deleteNodeContent(force); }
        free(d->pChild);
        for(i=0; i<d->nText; i++) free((void*)d->pText[i]);
        free(d->pText);
        for(i=0; i<d->nClear; i++) free((void*)d->pClear[i].lpszValue);
        free(d->pClear);
        for(i=0; i<d->nAttribute; i++)
        {
            free((void*)d->pAttribute[i].lpszName);
            if (d->pAttribute[i].lpszValue) free((void*)d->pAttribute[i].lpszValue);
        }
        free(d->pAttribute);
        free(d->pOrder);
        free((void*)d->lpszName);
        free(d);
        d=NULL;
    }
}

XMLNode XMLNode::addChild(XMLNode childNode, int pos)
{
    XMLNodeData *dc=childNode.d;
    if ((!dc)||(!d)) return childNode;
    if (dc->pParent) { if ((detachFromParent(dc)<=pos)&&(dc->pParent==d)) pos--; } else dc->ref_count++;
    dc->pParent=d;
//     int nc=d->nChild;
//     d->pChild=(XMLNode*)myRealloc(d->pChild,(nc+1),memoryIncrease,sizeof(XMLNode));
    d->pChild=(XMLNode*)addToOrder(0,&pos,d->nChild,d->pChild,sizeof(XMLNode),eNodeChild);
    d->pChild[pos].d=dc;
    d->nChild++;
    return childNode;
}

void XMLNode::deleteAttribute(int i)
{
    if ((!d)||(i<0)||(i>=d->nAttribute)) return;
    d->nAttribute--;
    XMLAttribute *p=d->pAttribute+i;
    free((void*)p->lpszName);
    if (p->lpszValue) free((void*)p->lpszValue);
    if (d->nAttribute) memmove(p,p+1,(d->nAttribute-i)*sizeof(XMLAttribute)); else { free(p); d->pAttribute=NULL; }
}

void XMLNode::deleteAttribute(XMLAttribute *a){ if (a) deleteAttribute(a->lpszName); }
void XMLNode::deleteAttribute(XMLCSTR lpszName)
{
    int j=0;
    getAttribute(lpszName,&j);
    if (j) deleteAttribute(j-1);
}

XMLAttribute *XMLNode::updateAttribute_WOSD(XMLCSTR lpszNewValue, XMLCSTR lpszNewName,int i)
{
    if (!d) return NULL;
    if (i>=d->nAttribute)
    {
        if (lpszNewName) return addAttribute_WOSD(lpszNewName,lpszNewValue);
        return NULL;
    }
    XMLAttribute *p=d->pAttribute+i;
    if (p->lpszValue&&p->lpszValue!=lpszNewValue) free((void*)p->lpszValue);
    p->lpszValue=lpszNewValue;
    if (lpszNewName&&p->lpszName!=lpszNewName) { free((void*)p->lpszName); p->lpszName=lpszNewName; };
    return p;
}

XMLAttribute *XMLNode::updateAttribute_WOSD(XMLAttribute *newAttribute, XMLAttribute *oldAttribute)
{
    if (oldAttribute) return updateAttribute_WOSD(newAttribute->lpszValue,newAttribute->lpszName,oldAttribute->lpszName);
    return addAttribute_WOSD(newAttribute->lpszName,newAttribute->lpszValue);
}

XMLAttribute *XMLNode::updateAttribute_WOSD(XMLCSTR lpszNewValue, XMLCSTR lpszNewName,XMLCSTR lpszOldName)
{
    int j=0;
    getAttribute(lpszOldName,&j);
    if (j) return updateAttribute_WOSD(lpszNewValue,lpszNewName,j-1);
    else
    {
        if (lpszNewName) return addAttribute_WOSD(lpszNewName,lpszNewValue);
        else             return addAttribute_WOSD(stringDup(lpszOldName),lpszNewValue);
    }
}

int XMLNode::indexText(XMLCSTR lpszValue) const
{
    if (!d) return -1;
    int i,l=d->nText;
    if (!lpszValue) { if (l) return 0; return -1; }
    XMLCSTR *p=d->pText;
    for (i=0; i<l; i++) if (lpszValue==p[i]) return i;
    return -1;
}

void XMLNode::deleteText(int i)
{
    if ((!d)||(i<0)||(i>=d->nText)) return;
    d->nText--;
    XMLCSTR *p=d->pText+i;
    free((void*)*p);
    if (d->nText) memmove(p,p+1,(d->nText-i)*sizeof(XMLCSTR)); else { free(p); d->pText=NULL; }
    removeOrderElement(d,eNodeText,i);
}

void XMLNode::deleteText(XMLCSTR lpszValue) { deleteText(indexText(lpszValue)); }

XMLCSTR XMLNode::updateText_WOSD(XMLCSTR lpszNewValue, int i)
{
    if (!d) return NULL;
    if (i>=d->nText) return addText_WOSD(lpszNewValue);
    XMLCSTR *p=d->pText+i;
    if (*p!=lpszNewValue) { free((void*)*p); *p=lpszNewValue; }
    return lpszNewValue;
}

XMLCSTR XMLNode::updateText_WOSD(XMLCSTR lpszNewValue, XMLCSTR lpszOldValue)
{
    if (!d) return NULL;
    int i=indexText(lpszOldValue);
    if (i>=0) return updateText_WOSD(lpszNewValue,i);
    return addText_WOSD(lpszNewValue);
}

void XMLNode::deleteClear(int i)
{
    if ((!d)||(i<0)||(i>=d->nClear)) return;
    d->nClear--;
    XMLClear *p=d->pClear+i;
    free((void*)p->lpszValue);
    if (d->nClear) memmove(p,p+1,(d->nClear-i)*sizeof(XMLClear)); else { free(p); d->pClear=NULL; }
    removeOrderElement(d,eNodeClear,i);
}

int XMLNode::indexClear(XMLCSTR lpszValue) const
{
    if (!d) return -1;
    int i,l=d->nClear;
    if (!lpszValue) { if (l) return 0; return -1; }
    XMLClear *p=d->pClear;
    for (i=0; i<l; i++) if (lpszValue==p[i].lpszValue) return i;
    return -1;
}

void XMLNode::deleteClear(XMLCSTR lpszValue) { deleteClear(indexClear(lpszValue)); }
void XMLNode::deleteClear(XMLClear *a) { if (a) deleteClear(a->lpszValue); }

XMLClear *XMLNode::updateClear_WOSD(XMLCSTR lpszNewContent, int i)
{
    if (!d) return NULL;
    if (i>=d->nClear)
    {
        return addClear_WOSD(XMLClearTags[0].lpszOpen,lpszNewContent,XMLClearTags[0].lpszClose);
    }
    XMLClear *p=d->pClear+i;
    if (lpszNewContent!=p->lpszValue) { free((void*)p->lpszValue); p->lpszValue=lpszNewContent; }
    return p;
}

XMLClear *XMLNode::updateClear_WOSD(XMLCSTR lpszNewValue, XMLCSTR lpszOldValue)
{
    if (!d) return NULL;
    int i=indexClear(lpszOldValue);
    if (i>=0) return updateClear_WOSD(lpszNewValue,i);
    return addClear_WOSD(lpszNewValue,XMLClearTags[0].lpszOpen,XMLClearTags[0].lpszClose);
}

XMLClear *XMLNode::updateClear_WOSD(XMLClear *newP,XMLClear *oldP)
{
    if (oldP) return updateClear_WOSD(newP->lpszValue,oldP->lpszValue);
    return NULL;
}

XMLNode& XMLNode::operator=( const XMLNode& A )
{
    // shallow copy
    if (this != &A)
    {
        deleteNodeContent();
        d=A.d;
        if (d) (d->ref_count) ++ ;
    }
    return *this;
}

XMLNode::XMLNode(const XMLNode &A)
{
    // shallow copy
    d=A.d;
    if (d) (d->ref_count)++ ;
}

int XMLNode::nChildNode(XMLCSTR name) const
{
    if (!d) return 0;
    int i,j=0,n=d->nChild;
    XMLNode *pc=d->pChild;
    for (i=0; i<n; i++)
    {
        if (_tcsicmp(pc->d->lpszName, name)==0) j++;
        pc++;
    }
    return j;
}

XMLNode XMLNode::getChildNode(XMLCSTR name, int *j) const
{
    if (!d) return emptyXMLNode;
    int i=0,n=d->nChild;
    if (j) i=*j;
    XMLNode *pc=d->pChild+i;
    for (; i<n; i++)
    {
        if (_tcsicmp(pc->d->lpszName, name)==0)
        {
            if (j) *j=i+1;
            return *pc;
        }
        pc++;
    }
    return emptyXMLNode;
}

XMLNode XMLNode::getChildNode(XMLCSTR name, int j) const
{
    if (!d) return emptyXMLNode;
    int i=0;
    while (j-->0) getChildNode(name,&i);
    return getChildNode(name,&i);
}

int XMLNode::positionOfText     (int i) const { if (i>=d->nText ) i=d->nText-1;  return findPosition(d,i,eNodeText ); }
int XMLNode::positionOfClear    (int i) const { if (i>=d->nClear) i=d->nClear-1; return findPosition(d,i,eNodeClear); }
int XMLNode::positionOfChildNode(int i) const { if (i>=d->nChild) i=d->nChild-1; return findPosition(d,i,eNodeChild); }
int XMLNode::positionOfText (XMLCSTR lpszValue) const { return positionOfText (indexText (lpszValue)); }
int XMLNode::positionOfClear(XMLCSTR lpszValue) const { return positionOfClear(indexClear(lpszValue)); }
int XMLNode::positionOfClear(XMLClear *a) const { if (a) return positionOfClear(a->lpszValue); return positionOfClear(); }
int XMLNode::positionOfChildNode(XMLNode x)  const
{
    if ((!d)||(!x.d)) return -1;
    XMLNodeData *dd=x.d;
    XMLNode *pc=d->pChild;
    int i=d->nChild;
    while (i--) if (pc[i].d==dd) return findPosition(d,i,eNodeChild);
    return -1;
}
int XMLNode::positionOfChildNode(XMLCSTR name, int count) const
{
    if (!name) return positionOfChildNode(count);
    int j=0;
    do { getChildNode(name,&j); if (j<0) return -1; } while (count--);
    return findPosition(d,j-1,eNodeChild);
}

XMLNode XMLNode::getChildNodeWithAttribute(XMLCSTR name,XMLCSTR attributeName,XMLCSTR attributeValue, int *k) const
{
     int i=0,j;
     if (k) i=*k;
     XMLNode x;
     XMLCSTR t;
     do
     {
         x=getChildNode(name,&i);
         if (!x.isEmpty())
         {
             if (attributeValue)
             {
                 j=0;
                 do
                 {
                     t=x.getAttribute(attributeName,&j);
                     if (t&&(_tcsicmp(attributeValue,t)==0)) { if (k) *k=i+1; return x; }
                 } while (t);
             } else
             {
                 if (x.isAttributeSet(attributeName)) { if (k) *k=i+1; return x; }
             }
         }
     } while (!x.isEmpty());
     return emptyXMLNode;
}

// Find an attribute on an node.
XMLCSTR XMLNode::getAttribute(XMLCSTR lpszAttrib, int *j) const
{
    if (!d) return NULL;
    int i=0,n=d->nAttribute;
    if (j) i=*j;
    XMLAttribute *pAttr=d->pAttribute+i;
    for (; i<n; i++)
    {
        if (_tcsicmp(pAttr->lpszName, lpszAttrib)==0)
        {
            if (j) *j=i+1;
            return pAttr->lpszValue;
        }
        pAttr++;
    }
    return NULL;
}

char XMLNode::isAttributeSet(XMLCSTR lpszAttrib) const
{
    if (!d) return FALSE;
    int i,n=d->nAttribute;
    XMLAttribute *pAttr=d->pAttribute;
    for (i=0; i<n; i++)
    {
        if (_tcsicmp(pAttr->lpszName, lpszAttrib)==0)
        {
            return TRUE;
        }
        pAttr++;
    }
    return FALSE;
}

XMLCSTR XMLNode::getAttribute(XMLCSTR name, int j) const
{
    if (!d) return NULL;
    int i=0;
    while (j-->0) getAttribute(name,&i);
    return getAttribute(name,&i);
}

XMLNodeContents XMLNode::enumContents(int i) const
{
    XMLNodeContents c;
    if (!d) { c.type=eNodeNULL; return c; }
    if (i<d->nAttribute)
    {
        c.type=eNodeAttribute;
        c.attrib=d->pAttribute[i];
        return c;
    }
    i-=d->nAttribute;
    c.type=(XMLElementType)(d->pOrder[i]&3);
    i=(d->pOrder[i])>>2;
    switch (c.type)
    {
    case eNodeChild:     c.child = d->pChild[i];      break;
    case eNodeText:      c.text  = d->pText[i];       break;
    case eNodeClear:     c.clear = d->pClear[i];      break;
    default: break;
    }
    return c;
}

XMLCSTR XMLNode::getName() const { if (!d) return NULL; return d->lpszName;   }
int XMLNode::nText()       const { if (!d) return 0;    return d->nText;      }
int XMLNode::nChildNode()  const { if (!d) return 0;    return d->nChild;     }
int XMLNode::nAttribute()  const { if (!d) return 0;    return d->nAttribute; }
int XMLNode::nClear()      const { if (!d) return 0;    return d->nClear;     }
int XMLNode::nElement()    const { if (!d) return 0;    return d->nAttribute+d->nChild+d->nText+d->nClear; }
XMLClear     XMLNode::getClear         (int i) const { if ((!d)||(i>=d->nClear    )) return emptyXMLClear;     return d->pClear[i];     }
XMLAttribute XMLNode::getAttribute     (int i) const { if ((!d)||(i>=d->nAttribute)) return emptyXMLAttribute; return d->pAttribute[i]; }
XMLCSTR      XMLNode::getAttributeName (int i) const { if ((!d)||(i>=d->nAttribute)) return NULL;              return d->pAttribute[i].lpszName;  }
XMLCSTR      XMLNode::getAttributeValue(int i) const { if ((!d)||(i>=d->nAttribute)) return NULL;              return d->pAttribute[i].lpszValue; }
XMLCSTR      XMLNode::getText          (int i) const { if ((!d)||(i>=d->nText     )) return NULL;              return d->pText[i];      }
XMLNode      XMLNode::getChildNode     (int i) const { if ((!d)||(i>=d->nChild    )) return emptyXMLNode;      return d->pChild[i];     }
XMLNode      XMLNode::getParentNode    (     ) const { if ((!d)||(!d->pParent     )) return emptyXMLNode;      return XMLNode(d->pParent); }
char         XMLNode::isDeclaration    (     ) const { if (!d) return 0;             return d->isDeclaration; }
char         XMLNode::isEmpty          (     ) const { return (d==NULL); }

XMLNode       XMLNode::addChild(XMLCSTR lpszName, char isDeclaration, int pos)
              { return addChild_priv(0,stringDup(lpszName),isDeclaration,pos); }
XMLNode       XMLNode::addChild_WOSD(XMLCSTR lpszName, char isDeclaration, int pos)
              { return addChild_priv(0,lpszName,isDeclaration,pos); }
XMLAttribute *XMLNode::addAttribute(XMLCSTR lpszName, XMLCSTR lpszValue)
              { return addAttribute_priv(0,stringDup(lpszName),stringDup(lpszValue)); }
XMLAttribute *XMLNode::addAttribute_WOSD(XMLCSTR lpszName, XMLCSTR lpszValuev)
              { return addAttribute_priv(0,lpszName,lpszValuev); }
XMLCSTR       XMLNode::addText(XMLCSTR lpszValue, int pos)
              { return addText_priv(0,stringDup(lpszValue),pos); }
XMLCSTR       XMLNode::addText_WOSD(XMLCSTR lpszValue, int pos)
              { return addText_priv(0,lpszValue,pos); }
XMLClear     *XMLNode::addClear(XMLCSTR lpszValue, XMLCSTR lpszOpen, XMLCSTR lpszClose, int pos)
              { return addClear_priv(0,stringDup(lpszValue),lpszOpen,lpszClose,pos); }
XMLClear     *XMLNode::addClear_WOSD(XMLCSTR lpszValue, XMLCSTR lpszOpen, XMLCSTR lpszClose, int pos)
              { return addClear_priv(0,lpszValue,lpszOpen,lpszClose,pos); }
XMLCSTR       XMLNode::updateName(XMLCSTR lpszName)
              { return updateName_WOSD(stringDup(lpszName)); }
XMLAttribute *XMLNode::updateAttribute(XMLAttribute *newAttribute, XMLAttribute *oldAttribute)
              { return updateAttribute_WOSD(stringDup(newAttribute->lpszValue),stringDup(newAttribute->lpszName),oldAttribute->lpszName); }
XMLAttribute *XMLNode::updateAttribute(XMLCSTR lpszNewValue, XMLCSTR lpszNewName,int i)
              { return updateAttribute_WOSD(stringDup(lpszNewValue),stringDup(lpszNewName),i); }
XMLAttribute *XMLNode::updateAttribute(XMLCSTR lpszNewValue, XMLCSTR lpszNewName,XMLCSTR lpszOldName)
              { return updateAttribute_WOSD(stringDup(lpszNewValue),stringDup(lpszNewName),lpszOldName); }
XMLCSTR       XMLNode::updateText(XMLCSTR lpszNewValue, int i)
              { return updateText_WOSD(stringDup(lpszNewValue),i); }
XMLCSTR       XMLNode::updateText(XMLCSTR lpszNewValue, XMLCSTR lpszOldValue)
              { return updateText_WOSD(stringDup(lpszNewValue),lpszOldValue); }
XMLClear     *XMLNode::updateClear(XMLCSTR lpszNewContent, int i)
              { return updateClear_WOSD(stringDup(lpszNewContent),i); }
XMLClear     *XMLNode::updateClear(XMLCSTR lpszNewValue, XMLCSTR lpszOldValue)
              { return updateClear_WOSD(stringDup(lpszNewValue),lpszOldValue); }
XMLClear     *XMLNode::updateClear(XMLClear *newP,XMLClear *oldP)
              { return updateClear_WOSD(stringDup(newP->lpszValue),oldP->lpszValue); }

void XMLNode::setGlobalOptions(char _guessUnicodeChars, char _strictUTF8Parsing, char _dropWhiteSpace)
{
    guessUnicodeChars=_guessUnicodeChars; dropWhiteSpace=_dropWhiteSpace; strictUTF8Parsing=_strictUTF8Parsing;
#ifndef _XMLUNICODE
    if (_strictUTF8Parsing) XML_ByteTable=XML_utf8ByteTable; else XML_ByteTable=XML_asciiByteTable;
#endif
}

char XMLNode::guessUTF8ParsingParameterValue(void *buf,int l, char useXMLEncodingAttribute)
{
#ifdef _XMLUNICODE
    return 0;
#else
    if (l<25) return 0;
    if (myIsTextUnicode(buf,l)) return 0;
    unsigned char *b=(unsigned char*)buf;
    if ((b[0]==0xef)&&(b[1]==0xbb)&&(b[2]==0xbf)) return 1;

    // Match utf-8 model ?
    int i=0;
    while (i<l)
        switch (XML_utf8ByteTable[b[i]])
        {
        case 4: i++; if ((i<l)&&(b[i]& 0xC0)!=0x80) return 0; // 10bbbbbb ?
        case 3: i++; if ((i<l)&&(b[i]& 0xC0)!=0x80) return 0; // 10bbbbbb ?
        case 2: i++; if ((i<l)&&(b[i]& 0xC0)!=0x80) return 0; // 10bbbbbb ?
        case 1: i++; break;
        case 0: i=l;
        }
    if (!useXMLEncodingAttribute) return 1;
    // if encoding is specified and different from utf-8 than it's non-utf8
    // otherwise it's utf-8
    char bb[201];
    l=mmin(l,200);
    memcpy(bb,buf,l); // copy buf into bb to be able to do "bb[l]=0"
    bb[l]=0;
    b=(unsigned char*)strstr(bb,"encoding");
    if (!b) return 1;
    b+=8; while XML_isSPACECHAR(*b) b++; if (*b!='=') return 1;
    b++;  while XML_isSPACECHAR(*b) b++; if ((*b!='\'')&&(*b!='"')) return 1;
    b++;  while XML_isSPACECHAR(*b) b++; if ((_strnicmp((char*)b,"utf-8",5)==0)||
                                             (_strnicmp((char*)b,"utf8",4)==0)) return 1;
    return 0;
#endif
}
#undef XML_isSPACECHAR

//////////////////////////////////////////////////////////
//      Here starts the base64 conversion functions.    //
//////////////////////////////////////////////////////////

static const char base64Fillchar = _T('='); // used to mark partial words at the end

// this lookup table defines the base64 encoding
XMLCSTR base64EncodeTable=_T("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/");

// Decode Table gives the index of any valid base64 character in the Base64 table]
// 96: '='  -   97: space char   -   98: illegal char   -   99: end of string
const unsigned char base64DecodeTable[] = {
    99,98,98,98,98,98,98,98,98,97,  97,98,98,97,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  //00 -29
    98,98,97,98,98,98,98,98,98,98,  98,98,98,62,98,98,98,63,52,53,  54,55,56,57,58,59,60,61,98,98,  //30 -59
    98,96,98,98,98, 0, 1, 2, 3, 4,   5, 6, 7, 8, 9,10,11,12,13,14,  15,16,17,18,19,20,21,22,23,24,  //60 -89
    25,98,98,98,98,98,98,26,27,28,  29,30,31,32,33,34,35,36,37,38,  39,40,41,42,43,44,45,46,47,48,  //90 -119
    49,50,51,98,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  //120 -149
    98,98,98,98,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  //150 -179
    98,98,98,98,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  //180 -209
    98,98,98,98,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  98,98,98,98,98,98,98,98,98,98,  //210 -239
    98,98,98,98,98,98,98,98,98,98,  98,98,98,98,98,98                                               //240 -255
};

XMLParserBase64Tool::~XMLParserBase64Tool(){ freeBuffer(); }

void XMLParserBase64Tool::freeBuffer(){ if (buf) free(buf); buf=NULL; buflen=0; }

int XMLParserBase64Tool::encodeLength(int inlen, char formatted)
{
    unsigned int i=((inlen-1)/3*4+4+1);
    if (formatted) i+=inlen/54;
    return i;
}

XMLSTR XMLParserBase64Tool::encode(unsigned char *inbuf, unsigned int inlen, char formatted)
{
    int i=encodeLength(inlen,formatted),k=17,eLen=inlen/3,j;
    alloc(i*sizeof(XMLCHAR));
    XMLSTR curr=(XMLSTR)buf;
    for(i=0;i<eLen;i++)
    {
        // Copy next three bytes into lower 24 bits of int, paying attention to sign.
        j=(inbuf[0]<<16)|(inbuf[1]<<8)|inbuf[2]; inbuf+=3;
        // Encode the int into four chars
        *(curr++)=base64EncodeTable[ j>>18      ];
        *(curr++)=base64EncodeTable[(j>>12)&0x3f];
        *(curr++)=base64EncodeTable[(j>> 6)&0x3f];
        *(curr++)=base64EncodeTable[(j    )&0x3f];
        if (formatted) { if (!k) { *(curr++)=_T('\n'); k=18; } k--; }
    }
    eLen=inlen-eLen*3; // 0 - 2.
    if (eLen==1)
    {
        *(curr++)=base64EncodeTable[ inbuf[0]>>2      ];
        *(curr++)=base64EncodeTable[(inbuf[0]<<4)&0x3F];
        *(curr++)=base64Fillchar;
        *(curr++)=base64Fillchar;
    } else if (eLen==2)
    {
        j=(inbuf[0]<<8)|inbuf[1];
        *(curr++)=base64EncodeTable[ j>>10      ];
        *(curr++)=base64EncodeTable[(j>> 4)&0x3f];
        *(curr++)=base64EncodeTable[(j<< 2)&0x3f];
        *(curr++)=base64Fillchar;
    }
    *(curr++)=0;
    return (XMLSTR)buf;
}

unsigned int XMLParserBase64Tool::decodeSize(XMLCSTR data,XMLError *xe)
{
     if (xe) *xe=eXMLErrorNone;
    int size=0;
    unsigned char c;
    //skip any extra characters (e.g. newlines or spaces)
    while (*data)
    {
#ifdef _XMLUNICODE
        if (*data>255) { if (xe) *xe=eXMLErrorBase64DecodeIllegalCharacter; return 0; }
#endif
        c=base64DecodeTable[(unsigned char)(*data)];
        if (c<97) size++;
        else if (c==98) { if (xe) *xe=eXMLErrorBase64DecodeIllegalCharacter; return 0; }
        data++;
    }
    if (xe&&(size%4!=0)) *xe=eXMLErrorBase64DataSizeIsNotMultipleOf4;
    if (size==0) return 0;
    do { data--; size--; } while(*data==base64Fillchar); size++;
    return (unsigned int)((size*3)/4);
}

unsigned char XMLParserBase64Tool::decode(XMLCSTR data, unsigned char *buf, int len, XMLError *xe)
{
    if (xe) *xe=eXMLErrorNone;
    int i=0,p=0;
    unsigned char d,c;
    for(;;)
    {

#ifdef _XMLUNICODE
#define BASE64DECODE_READ_NEXT_CHAR(c)                                              \
        do {                                                                        \
            if (data[i]>255){ c=98; break; }                                        \
            c=base64DecodeTable[(unsigned char)data[i++]];                       \
        }while (c==97);                                                             \
        if(c==98){ if(xe)*xe=eXMLErrorBase64DecodeIllegalCharacter; return 0; }
#else
#define BASE64DECODE_READ_NEXT_CHAR(c)                                           \
        do { c=base64DecodeTable[(unsigned char)data[i++]]; }while (c==97);   \
        if(c==98){ if(xe)*xe=eXMLErrorBase64DecodeIllegalCharacter; return 0; }
#endif

        BASE64DECODE_READ_NEXT_CHAR(c)
        if (c==99) { return 2; }
        if (c==96)
        {
            if (p==(int)len) return 2;
            if (xe) *xe=eXMLErrorBase64DecodeTruncatedData;
            return 1;
        }

        BASE64DECODE_READ_NEXT_CHAR(d)
        if ((d==99)||(d==96)) { if (xe) *xe=eXMLErrorBase64DecodeTruncatedData;  return 1; }
        if (p==(int)len) {      if (xe) *xe=eXMLErrorBase64DecodeBufferTooSmall; return 0; }
        buf[p++]=(c<<2)|((d>>4)&0x3);

        BASE64DECODE_READ_NEXT_CHAR(c)
        if (c==99) { if (xe) *xe=eXMLErrorBase64DecodeTruncatedData;  return 1; }
        if (p==(int)len)
        {
            if (c==96) return 2;
            if (xe) *xe=eXMLErrorBase64DecodeBufferTooSmall;
            return 0;
        }
        if (c==96) { if (xe) *xe=eXMLErrorBase64DecodeTruncatedData;  return 1; }
        buf[p++]=((d<<4)&0xf0)|((c>>2)&0xf);

        BASE64DECODE_READ_NEXT_CHAR(d)
        if (d==99 ) { if (xe) *xe=eXMLErrorBase64DecodeTruncatedData;  return 1; }
        if (p==(int)len)
        {
            if (d==96) return 2;
            if (xe) *xe=eXMLErrorBase64DecodeBufferTooSmall;
            return 0;
        }
        if (d==96) { if (xe) *xe=eXMLErrorBase64DecodeTruncatedData;  return 1; }
        buf[p++]=((c<<6)&0xc0)|d;
    }
}
#undef BASE64DECODE_READ_NEXT_CHAR

void XMLParserBase64Tool::alloc(int newsize)
{
    if ((!buf)&&(newsize)) { buf=malloc(newsize); buflen=newsize; return; }
    if (newsize>buflen) { buf=realloc(buf,newsize); buflen=newsize; }
}

unsigned char *XMLParserBase64Tool::decode(XMLCSTR data, int *outlen, XMLError *xe)
{
    if (xe) *xe=eXMLErrorNone;
    unsigned int len=decodeSize(data,xe);
    if (outlen) *outlen=len;
    if (!len) return NULL;
    alloc(len+1);
    if(!decode(data,(unsigned char*)buf,len,xe)){ return NULL; }
    return (unsigned char*)buf;
}

