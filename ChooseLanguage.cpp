// ChooseLanguage.cpp : implementation file
//

#include "stdafx.h"
#include "OCRProject.h"
#include "ChooseLanguage.h"
#include "afxdialogex.h"


// CChooseLanguage dialog

IMPLEMENT_DYNAMIC(CChooseLanguage, CDialogEx)

CChooseLanguage::CChooseLanguage(CWnd* pParent /*=NULL*/)
	: CDialogEx(IDD_DIALOG_CHOOSE_LANGUAGE, pParent)
	, m_language(0)
{

}

CChooseLanguage::~CChooseLanguage()
{
}

void CChooseLanguage::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Radio(pDX, IDC_RADIO_CHINESE, m_language);
}


BEGIN_MESSAGE_MAP(CChooseLanguage, CDialogEx)
//	ON_BN_CLICKED(IDC_RADIO_CHINESE, &CChooseLanguage::OnBnClickedOk)
	ON_BN_CLICKED(IDC_RADIO_ENGLISH, &CChooseLanguage::OnBnClickedRadioChinese)
	ON_BN_CLICKED(IDC_RADIO_CHINESE, &CChooseLanguage::OnBnClickedRadioChinese)
	ON_BN_CLICKED(IDOK, &CChooseLanguage::OnBnClickedOk)
END_MESSAGE_MAP()


// CChooseLanguage message handlers


//void CChooseLanguage::OnBnClickedOk()
//{
//	UpdateData(true);
//	//while (true)
//	//{
//	//	switch (m_language)
//	//	{
//	//	//case 0:MessageBox(_T("ÖÐÎÄ")); break;
//	//	case 1:MessageBox(_T("Ó¢ÎÄ")); break;
//	//	}
//	//}
//	//CDialogEx::OnOK();
//}


void CChooseLanguage::OnBnClickedRadioChinese()
{
	// TODO: Add your control notification handler code here
}


void CChooseLanguage::OnBnClickedOk()
{
	CDialogEx::OnOK();
}
