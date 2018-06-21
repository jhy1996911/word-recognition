#pragma once
#include "afxwin.h"


// CChooseLanguage dialog

class CChooseLanguage : public CDialogEx
{
	DECLARE_DYNAMIC(CChooseLanguage)

public:
	CChooseLanguage(CWnd* pParent = NULL);   // standard constructor
	virtual ~CChooseLanguage();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_DIALOG_CHOOSE_LANGUAGE };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	int m_language=0;
//	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedRadioChinese();
	afx_msg void OnBnClickedOk();
	CEdit m_tessdata;
};
