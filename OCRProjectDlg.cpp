
// OCRProjectDlg.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "OCRProject.h"
#include "OCRProjectDlg.h"
#include "afxdialogex.h"
#include<iostream>    
#include <opencv2/core/core.hpp>    
#include <opencv2/highgui/highgui.hpp> 
#include <cstring> 
#include <opencv2/imgproc/imgproc.hpp>  
#include <baseapi.h>  
#include <iostream> 
#include <cstdlib>
#include "PreDeal.h"
#include "FindLocation.h"
#include "ChooseLanguage.h"
using namespace cv;
using namespace std;
#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// COCRProjectDlg �Ի���

COCRProjectDlg::COCRProjectDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(IDD_OCRPROJECT_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void COCRProjectDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(COCRProjectDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON1, &COCRProjectDlg::OnBnClickedButton1)
	ON_COMMAND(ID_32771, &COCRProjectDlg::open_img)
	ON_BN_CLICKED(IDC_BUTTON2, &COCRProjectDlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_BUTTON3, &COCRProjectDlg::OnBnClickedButton3)
	ON_COMMAND(ID_32775, &COCRProjectDlg::Threshold)
	//ON_COMMAND(ID_32776, &COCRProjectDlg::blur)
	ON_COMMAND(ID_32780, &COCRProjectDlg::OnWriteText)
	ON_COMMAND(ID_32781, &COCRProjectDlg::OnBlackText)
	ON_COMMAND(ID_MSER_32783, &COCRProjectDlg::OnMserWhite)
	ON_COMMAND(ID_MSER_32784, &COCRProjectDlg::OnMserBlack)
	ON_COMMAND(ID_32776, &COCRProjectDlg::On32776)
	ON_COMMAND(ID_32785, &COCRProjectDlg::On32785)
	ON_COMMAND(ID_32786, &COCRProjectDlg::On32786)
	ON_COMMAND(ID_32787, &COCRProjectDlg::On32787)
	ON_COMMAND(ID_32772, &COCRProjectDlg::OnSavePic)
END_MESSAGE_MAP()


// COCRProjectDlg ��Ϣ�������

BOOL COCRProjectDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// ���ô˶Ի����ͼ�ꡣ  ��Ӧ�ó��������ڲ��ǶԻ���ʱ����ܽ��Զ�
	//  ִ�д˲���
	SetIcon(m_hIcon, TRUE);			// ���ô�ͼ��
	SetIcon(m_hIcon, FALSE);		// ����Сͼ��

	// TODO: �ڴ���Ӷ���ĳ�ʼ������

	return TRUE;  // ���ǽ��������õ��ؼ������򷵻� TRUE
}

// �����Ի��������С����ť������Ҫ����Ĵ���
//  �����Ƹ�ͼ�ꡣ  ����ʹ���ĵ�/��ͼģ�͵� MFC Ӧ�ó���
//  �⽫�ɿ���Զ���ɡ�

void COCRProjectDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // ���ڻ��Ƶ��豸������

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// ʹͼ���ڹ����������о���
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// ����ͼ��
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//���û��϶���С������ʱϵͳ���ô˺���ȡ�ù��
//��ʾ��
HCURSOR COCRProjectDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}


char *out;
Mat im;
vector<Rect>part;
bool FindDeal=false;
Mat resPart[100];
void COCRProjectDlg::open_img()
{
	FindDeal = false;
	CString filter;
	CFileDialog dlg(TRUE, NULL, NULL, OFN_HIDEREADONLY);
	CString imgPath;
	Mat sh;
	if (dlg.DoModal() == IDOK)
	{
		imgPath = dlg.GetPathName();
		UpdateData(FALSE);
	}
	USES_CONVERSION;
	// TODO: �ڴ���ӿؼ�֪ͨ����������
	string srcPath(W2A(imgPath.GetBuffer()));
	im = imread(srcPath);
	if (im.empty())
	{
		//std::cout << "Cannot open source image!" << std::endl;
		//MessageBox(_T("Cannot open source image!"));
		MessageBox(_T("û������ͼƬ"));
		return;
	}
	if(im.cols>1000&&im.rows>1000)
	resize(im,im, cv::Size(0, 0), 1.0/(1.0*im.cols/1000), 1.0/(1.0*im.rows/1000), cv::INTER_LINEAR);
	else if (im.cols < 500 && im.rows < 500)
	{
		resize(im, im, cv::Size(0, 0), 2, 2, cv::INTER_LINEAR);
	}
	cvNamedWindow("����");
	imshow("����", im);
	//cv::cvtColor(im, im, CV_BGR2GRAY);
	//equalizeHist(im, im);
	waitKey(0);
}

CString UTF82WCS(const char* szU8)
{
	//Ԥת�����õ�����ռ�Ĵ�С;
	int wcsLen = ::MultiByteToWideChar(CP_UTF8, NULL, szU8, strlen(szU8), NULL, 0);

	//����ռ�Ҫ��'\0'�����ռ䣬MultiByteToWideChar�����'\0'�ռ�
	wchar_t* wszString = new wchar_t[wcsLen + 1];

	//ת��
	::MultiByteToWideChar(CP_UTF8, NULL, szU8, strlen(szU8), wszString, wcsLen);

	//������'\0'
	wszString[wcsLen] = '\0';

	CString unicodeString(wszString);

	delete[] wszString;
	wszString = NULL;

	return unicodeString;
}

void arrayCopy(char *from, char * to)
{   //��Ҫ���׸ı��βε�ֵ�����ø���ָ����������βνӹ�����  
	char *myfrom = from;
	char *myto = to;
	//if (from == NULL || to == NULL)//if�����ǿ����׳��  
	//{
	//	return -1;
	//}
	while (myfrom != NULL)
	{
		*myto = *myfrom;
		myfrom++;
		myto++;
	}
}
void COCRProjectDlg::OnBnClickedButton1()
{
	CChooseLanguage chooselanguage;
	CString lin;
	chooselanguage.DoModal();
	char *chice=NULL;
	if (chooselanguage.m_language == 0)
		chice = "chi_sim";
	else if (chooselanguage.m_language == 1)
		chice = "eng";
	else if (chooselanguage.m_language == 2)
	{
		chooselanguage.m_tessdata.GetWindowText(lin);
		chice = (LPSTR)(LPCTSTR)lin;
	}
	if (im.empty())
	{
		//std::cout << "Cannot open source image!" << std::endl;
		//MessageBox(_T("Cannot open source image!"));
		MessageBox(_T("û������ͼƬ"));
		return;
	}	
	if (im.channels()!=1)
	{
		cv::cvtColor(im, im, CV_BGR2GRAY);
	}
	// ...other im.getMat()age pre-processing here...
	
	// Pass it to Tesseract API  

	CString str;
	if (FindDeal == false)
	{
		tesseract::TessBaseAPI tess;
		//tess.Init(NULL, "chi_sim", tesseract::OEM_DEFAULT);
		tess.Init(NULL, chice, tesseract::OEM_DEFAULT);
		//tess.SetPageSegMode(tesseract::PSM_SINGLE_BLOCK);
		tess.SetImage((uchar*)im.data, im.cols, im.rows, 1, im.cols);
		// Get the text  
		out = tess.GetUTF8Text();
		CEdit* pBoxOne;
		pBoxOne = (CEdit*)GetDlgItem(IDC_EDIT1);
		CString str1 = UTF82WCS(out);
		//��TCHAR����ת��ΪCString
		//ɾ��������
		for (int i = 0; i < str1.GetLength(); ++i)
		{
			if (i % 55 == 0 && i != 0)
				str += "\r\n";
			str += str1[i];
		}
		pBoxOne->SetWindowTextW(str);
	}
	else
	{
		for (int i = 0; i < part.size(); ++i)
		{
			tesseract::TessBaseAPI tess;
			//tess.Init(NULL, "chi_sim", tesseract::OEM_DEFAULT);
			tess.Init(NULL, chice, tesseract::OEM_DEFAULT);
			tess.SetPageSegMode(tesseract::PSM_SINGLE_BLOCK);
			//imshow("1", resPart[i]); 
			//waitKey(0);
			char lin[100];
			sprintf_s(lin,"deal/%d.jpg", i);
			imwrite(lin, resPart[i]);
			Mat a = imread(lin);
			imshow("1", a);
			if (a.channels() != 1)
			{
				cv::cvtColor(a, a, CV_BGR2GRAY);
			}
			waitKey(0);
			//tess.SetImage((uchar*)resPart[i].data, resPart[i].cols, resPart[i].rows, 1, resPart[i].cols);
			tess.SetImage((uchar*)a.data, a.cols, a.rows, 1, a.cols);
			// Get the text 
			//imshow("1", resPart[i]); 
			//waitKey(0);
			//if(tess.GetUTF8Text()!=NULL)
			out = tess.GetUTF8Text();	

			CString str1 = UTF82WCS(out);
				for (int j = 0; j < str1.GetLength(); ++j)
				{
					if (j % 55 == 0 && j != 0)
						str += "\r\n";
					str += str1[j];
				}
				str += "\r\n";
			
		}
		CEdit* pBoxOne;
		pBoxOne = (CEdit*)GetDlgItem(IDC_EDIT1); 
		pBoxOne->SetWindowTextW(str);

		
	}

}


void COCRProjectDlg::OnBnClickedButton2()
{
	CFile test;
	CFileException e;
	CFileDialog dlg(TRUE, NULL, NULL, OFN_HIDEREADONLY);
	CString savePath;
	if (dlg.DoModal() == IDOK)
	{
		savePath = dlg.GetPathName();
		UpdateData(FALSE);
	}
	TCHAR* pszFileName = (LPTSTR)(LPCTSTR)savePath;
	if (!test.Open(pszFileName, CFile::modeCreate | CFile::modeNoTruncate | CFile::modeReadWrite, &e))//��������test.txt�ļ�
	{
		MessageBox(_T("File could not be opened"));
		return;
	}
	CString m_str;
	GetDlgItemText(IDC_EDIT1, m_str); 
	test.Write(out,strlen(out));
	test.Close();//�ر��ļ�
}

void COCRProjectDlg::OnBnClickedButton3()
{
	if (im.empty())
	{
		std::cout << "Cannot open source image!" << std::endl;
		//MessageBox(_T("Cannot open source image!"));
		MessageBox(_T("û������ͼƬ"));
		return;
	}
	Mat showIm;
	showIm = im.clone();
	//if(showIm.cols>1000&& showIm.rows>1000)
	//resize(showIm, showIm, cv::Size(0, 0), 1.0/(1.0*showIm.cols/1000), 1.0/(1.0*showIm.rows/1000), cv::INTER_LINEAR);
	//else if (showIm.cols < 500 && showIm.rows < 500)
	//{
	//	resize(showIm, showIm, cv::Size(0, 0), 2, 2, cv::INTER_LINEAR);
	//}
	imshow("��ǰͼƬЧ��", showIm);
}
PreDeal deal;

int g_nMedianBlurValue = 5;  //��ֵ�˲�����ֵ




//��Ȼ�������ֶ�λ


int co = -1;



//����
void COCRProjectDlg::OnWriteText()
{
	co = 0;
	findlocation.mainTextDetection(im,co);
}

//����
void COCRProjectDlg::OnBlackText()
{
	co = 1;
	findlocation.mainTextDetection(im,co);
}



Mat reverse(Mat src)

{

	Mat dst = src<100;

	return dst;
}


void COCRProjectDlg::OnMserWhite()
{
	co = 0;
	findwordarea.FindArea(im,part);
	for (int i = 0; i < part.size(); ++i)
	{
		//resPart[i] = reverse(resPart[i]);
		resPart[i] = im(part[i]);
		cvtColor(resPart[i], resPart[i], CV_BGR2GRAY);
		//Deal.adaptiveThreshold(resPart[i], resPart[i], resPart[i].rows, resPart[i].cols);
		threshold(resPart[i], resPart[i], 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
		resPart[i] = reverse(resPart[i]);
	}
	FindDeal = true;
}





void COCRProjectDlg::OnMserBlack()
{
	co = 0;
	findwordarea.FindArea(im, part);
	for (int i = 0; i < part.size(); ++i)
	{
		resPart[i] = im(part[i]);
		if (resPart[i].channels() != 1)
		{
			cvtColor(resPart[i], resPart[i], CV_BGR2GRAY);
		}
		//Deal.adaptiveThreshold(resPart[i], resPart[i], resPart[i].rows, resPart[i].cols);
		threshold(resPart[i], resPart[i], 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
		//adaptiveThreshold(resPart[i], resPart[i], 255, CV_ADAPTIVE_THRESH_GAUSSIAN_C, CV_THRESH_BINARY, 9, 5);
	}
	FindDeal = true;
}


void COCRProjectDlg::On32776()
{
	GaussianBlur(im, im,Size(3,3),0,0);
	imshow("�˲�",im);
	waitKey(0);
}

int threshold_color = 0;

void COCRProjectDlg::On32785()
{
	cvtColor(im, im, COLOR_BGR2GRAY);
	imshow("�ҶȻ�ͼƬ",im);
	waitKey(0);
}

void COCRProjectDlg::Threshold()
{
	if (im.empty())
	{
		std::cout << "Cannot open source image!" << std::endl;
		//MessageBox(_T("Cannot open source image!"));
		MessageBox(_T("û������ͼƬ"));
		return;
	}
	if (2 == threshold_color)
	{
		im=reverse(im);
	}
	int width, height;
	if (im.cols < im.rows)
	{
		Deal.adaptiveThreshold(im, im, im.cols, im.rows);
		imshow("2", im);
		waitKey(0);
	}
	else
	{
		Deal.adaptiveThreshold(im, im, im.rows, im.cols);
	}
}


void COCRProjectDlg::On32786()
{
	threshold_color = 1;
	Threshold();
}


void COCRProjectDlg::On32787()
{
	threshold_color = 2;
	Threshold();
}

CString SelFilePath()
{
	TCHAR           szFolderPath[MAX_PATH] = { 0 };
	CString         strFolderPath = TEXT("");

	BROWSEINFO      sInfo;
	::ZeroMemory(&sInfo, sizeof(BROWSEINFO));
	sInfo.pidlRoot = 0;
	sInfo.lpszTitle = _T("��ѡ�������洢·��");
	sInfo.ulFlags = BIF_RETURNONLYFSDIRS | BIF_EDITBOX | BIF_DONTGOBELOWDOMAIN;
	sInfo.lpfn = NULL;

	// ��ʾ�ļ���ѡ��Ի���  
	LPITEMIDLIST lpidlBrowse = ::SHBrowseForFolder(&sInfo);
	if (lpidlBrowse != NULL)
	{
		// ȡ���ļ�����  
		if (::SHGetPathFromIDList(lpidlBrowse, szFolderPath))
		{
			strFolderPath = szFolderPath;
		}
	}
	if (lpidlBrowse != NULL)
	{
		::CoTaskMemFree(lpidlBrowse);
	}

	return strFolderPath;
}
void COCRProjectDlg::OnSavePic()
{
	// TODO: Add your command handler code here
	CString path=SelFilePath();
	String PicPath = CStringA(path);
	PicPath += "\\deal.jpg";
	imwrite(PicPath, im);

}
