
// OCRProject.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// COCRProjectApp: 
// �йش����ʵ�֣������ OCRProject.cpp
//

class COCRProjectApp : public CWinApp
{
public:
	COCRProjectApp();

// ��д
public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern COCRProjectApp theApp;
