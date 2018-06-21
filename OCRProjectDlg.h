
// OCRProjectDlg.h : 头文件
//

#pragma once
#include "Predeal.h"
#include "FindLocation.h"
#include "FindWordArea.h"
// COCRProjectDlg 对话框
class COCRProjectDlg : public CDialogEx
{
// 构造
public:
	COCRProjectDlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_OCRPROJECT_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButton1();
	afx_msg void open_img();
	afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedButton3();
	afx_msg void Threshold();
	afx_msg void OnWriteText();
	afx_msg void OnBlackText();
	afx_msg PreDeal Deal;
	afx_msg FindLocation findlocation;
	afx_msg FindWordArea findwordarea;
	afx_msg void OnMserWhite();
	afx_msg void OnMserBlack();
	afx_msg void On32776();
	afx_msg void On32785();
	afx_msg void On32786();
	afx_msg void On32787();
	afx_msg void OnSavePic();
};
