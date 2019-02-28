#include "stdafx.h"
#include <windows.h>
#include <d2d1.h>
#pragma comment(lib, "d2d1")

#include "basewin.h"
#include "matrix.h"


template <class T> void SafeRelease(T **ppT)
{
	if (*ppT)
	{
		(*ppT)->Release();
		*ppT = NULL;
	}
}


class MainWindow : public BaseWindow<MainWindow>
{
	const int regions = 20;
	const int correction = 10;

	ID2D1Factory            *pFactory;
	ID2D1HwndRenderTarget   *pRenderTarget;
	ID2D1SolidColorBrush    *pBrush;
	ID2D1SolidColorBrush 	*pathBrush;
	std::vector<ID2D1SolidColorBrush*> backBrushs;


	Matrix M;

	float xSize;
	float ySize;
	
	float xMid;
	float yMid;

	void    CalculateLayout();


	HRESULT CreateGraphicsResources();
	void    DiscardGraphicsResources();
	
	void    OnPaint();
	void    Resize();
	void	BackDrow();
	void	pathDrow();
	// перход от Графических координат к Декартовым
	// x = (Gx - xMid)/SIZECOEF <=> Gx = SIZECOEF*x + xMid
	// y = (yMid - Gy)/SIZECOEF <=> Gy = yMid - SIZECOEF*y

	const int SIZECOEF = 45;
	double Dx2Gx(double x) { return  SIZECOEF * x + xMid; };
	double Dy2Gy(double y) { return  yMid - SIZECOEF * y; };

	double Gx2Dx(double Gx) { return (Gx - xMid) / SIZECOEF; };
	double Gy2Dy(double Gy) { return (yMid - Gy) / SIZECOEF; };


public:

	MainWindow() : pFactory(NULL), pRenderTarget(NULL), pBrush(NULL)
	{
	}

	PCWSTR  ClassName() const { return L"Circle Window Class"; }
	LRESULT HandleMessage(UINT uMsg, WPARAM wParam, LPARAM lParam);
};


// Recalculate drawing layout when the size of the window changes.

void MainWindow::CalculateLayout()
{
	if (pRenderTarget != NULL)
	{
		D2D1_SIZE_F size = pRenderTarget->GetSize();

		xSize = size.width;
		ySize = size.height;

		xMid = xSize / 2 + 1;
		yMid = ySize / 2 + 1;
	}
}

HRESULT MainWindow::CreateGraphicsResources()
{
	HRESULT hr = S_OK;
	if (pRenderTarget == NULL)
	{
		RECT rc;
		GetClientRect(m_hwnd, &rc);

		D2D1_SIZE_U size = D2D1::SizeU(rc.right, rc.bottom);

		hr = pFactory->CreateHwndRenderTarget(
			D2D1::RenderTargetProperties(),
			D2D1::HwndRenderTargetProperties(m_hwnd, size),
			&pRenderTarget);

		if (SUCCEEDED(hr))
		{
			const D2D1_COLOR_F color = D2D1::ColorF(D2D1::ColorF::White);
			hr = pRenderTarget->CreateSolidColorBrush(color, &pBrush); /////////////////////////////////
			pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(D2D1::ColorF::Blue), &pathBrush);


			int frames = regions + correction;
			backBrushs.clear();
			backBrushs.resize(frames);

	
			for (int i = 0; i < frames; i++) {
				hr = pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(double(i) / frames,	double(i) / frames, double(i) / frames), &backBrushs[i]);
			}


			if (SUCCEEDED(hr))
			{
				CalculateLayout();
			}
		}
	}
	return hr;
}

void MainWindow::DiscardGraphicsResources()
{
	SafeRelease(&pRenderTarget);
	SafeRelease(&pBrush);
}


void MainWindow::OnPaint()					// Процесс отрисовки:
{
	HRESULT hr = CreateGraphicsResources();
	if (SUCCEEDED(hr))
	{
		PAINTSTRUCT ps;
		BeginPaint(m_hwnd, &ps);

		pRenderTarget->BeginDraw();		// Начало отрисовки
		BackDrow();						// Закраска фона

		// отрисовка заданных функций

#ifdef ONE
		D2D1_ELLIPSE   cyrcle1 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(0), Dy2Gy(0)), 1 * SIZECOEF, 1 * SIZECOEF);
		D2D1_ELLIPSE   cyrcle2 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(3), Dy2Gy(2)),  2 * SIZECOEF, 2 * SIZECOEF);

		pRenderTarget->DrawEllipse(cyrcle1, pBrush);
		pRenderTarget->DrawEllipse(cyrcle2, pBrush);
#endif

#ifdef TWO
		D2D1_ELLIPSE   cyrcle1 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(-2), Dy2Gy(-2)), 2 * SIZECOEF, 2 * SIZECOEF);
		D2D1_ELLIPSE   cyrcle2 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(1), Dy2Gy(2)), 3 * SIZECOEF, 3 * SIZECOEF);

		pRenderTarget->DrawEllipse(cyrcle1, pBrush);
		pRenderTarget->DrawEllipse(cyrcle2, pBrush);
#endif

#ifdef THREE
		D2D1_ELLIPSE   cyrcle1 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(-2), Dy2Gy(-2)), 3 * SIZECOEF, 3 * SIZECOEF);
		D2D1_ELLIPSE   cyrcle2 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(1), Dy2Gy(2)), 3 * SIZECOEF, 3 * SIZECOEF);

		pRenderTarget->DrawEllipse(cyrcle1, pBrush);
		pRenderTarget->DrawEllipse(cyrcle2, pBrush);
#endif


#ifdef  FOUR 
		D2D1_ELLIPSE   cyrcle1 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(0), Dy2Gy(0)), 1 * SIZECOEF, 1 * SIZECOEF);
		D2D1_ELLIPSE   cyrcle2 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(3), Dy2Gy(2)), 2 * SIZECOEF, 2 * SIZECOEF);
		D2D1_ELLIPSE   cyrcle3 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(3), Dy2Gy(-3)), 1.5 * SIZECOEF, 1.5 * SIZECOEF);

		pRenderTarget->DrawEllipse(cyrcle1, pBrush);
		pRenderTarget->DrawEllipse(cyrcle2, pBrush);
		pRenderTarget->DrawEllipse(cyrcle3, pBrush);
#endif

#ifdef  FIVE
		D2D1_ELLIPSE   cyrcle1 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(0), Dy2Gy(0)), 1 * SIZECOEF, 1 * SIZECOEF);
		D2D1_ELLIPSE   cyrcle2 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(3), Dy2Gy(2)), 2 * SIZECOEF, 2 * SIZECOEF);
		D2D1_ELLIPSE   cyrcle3 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(3), Dy2Gy(-3)), 1.5 * SIZECOEF, 1.5 * SIZECOEF);

		pRenderTarget->DrawEllipse(cyrcle1, pBrush);
		pRenderTarget->DrawEllipse(cyrcle2, pBrush);
		pRenderTarget->DrawEllipse(cyrcle3, pBrush);
#endif


#ifdef  SIX
		D2D1_ELLIPSE   cyrcle1 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(0), Dy2Gy(0)), 2 * SIZECOEF, 2 * SIZECOEF);
		D2D1_ELLIPSE   cyrcle2 = D2D1::Ellipse(D2D1::Point2F(Dx2Gx(4), Dy2Gy(0)), 2 * SIZECOEF, 2 * SIZECOEF);
		pRenderTarget->DrawEllipse(cyrcle1, pBrush);
		pRenderTarget->DrawEllipse(cyrcle2, pBrush);
#endif

		pathDrow();		// отрисовка пути решения

		////////////////////
		hr = pRenderTarget->EndDraw();	// Конец отрисовки
		if (FAILED(hr) || hr == D2DERR_RECREATE_TARGET)
		{
			DiscardGraphicsResources();
		}
		EndPaint(m_hwnd, &ps);
	}
}

void MainWindow::Resize()
{
	if (pRenderTarget != NULL)
	{
		RECT rc;
		GetClientRect(m_hwnd, &rc);

		D2D1_SIZE_U size = D2D1::SizeU(rc.right, rc.bottom);

		pRenderTarget->Resize(size);
		CalculateLayout();
		InvalidateRect(m_hwnd, NULL, FALSE);
	}
}

void MainWindow::BackDrow() {
	std::vector<double> regionsFrames(regions + correction);

	std::vector<double> f(M.getNumEquals());
	double maxNorm = 0;

	std::vector<std::vector<double> > backNorm(ySize);
	for (int i = 0; i < backNorm.size(); i++) backNorm[i].resize(xSize);
	///////////////////////////////////////////////////////////

	for (int y = 0; y < ySize; y++) {
		for (int x = 0; x < xSize; x++) {
			backNorm[y][x] = M.normInPoint(std::vector<double>{(x - xMid) / SIZECOEF, (yMid - y) / SIZECOEF }, f);
			if (backNorm[y][x] > maxNorm) maxNorm = backNorm[y][x];
		}
	}

	for (int i = 0; i < correction; i++) {
		regionsFrames[i] = double(i + 1) / (regions * correction) * maxNorm;
	}
	for (int i = correction - 1; i < regions + correction; i++) {
		regionsFrames[i] = double(i - correction + 2) / (regions)* maxNorm;
	}

	for (int y = 0; y < ySize; y++) {
		for (int x = 0; x < xSize; x++) {

			for (int region = 0; region < regions + correction; region++) {
				if (backNorm[y][x] < regionsFrames[region]) {
					D2D1_RECT_F rectangle = D2D1::RectF(x, y, x, y);
					pRenderTarget->DrawRectangle(&rectangle, backBrushs[region]); //
					break;	// loop exit
				}
			}
		}
	}

	// координатные линии
	pRenderTarget->DrawLine(D2D1::Point2F(0.0f, ySize / 2), D2D1::Point2F(xSize, ySize / 2), pBrush, 1.0f);
	pRenderTarget->DrawLine(D2D1::Point2F(xSize / 2, 0.0f), D2D1::Point2F(xSize / 2, ySize), pBrush, 1.0f);
}

void MainWindow::pathDrow() {
	M.solve();

	for (int i = 1; i < M.xPath.size(); i++) {
		pRenderTarget->DrawLine(D2D1::Point2F(SIZECOEF * M.xPath[i - 1][0] + xMid, yMid - SIZECOEF * M.xPath[i - 1][1]),
			D2D1::Point2F(SIZECOEF * M.xPath[i][0] + xMid, yMid - SIZECOEF * M.xPath[i][1]),
			pathBrush, 2.5f);
	}
}
int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE, PWSTR, int nCmdShow)
{
	MainWindow win;

	if (!win.Create(L"Circle", WS_OVERLAPPEDWINDOW))
	{
		return 0;
	}

	ShowWindow(win.Window(), nCmdShow);

	// Run the message loop.

	MSG msg = { };
	while (GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}

	return 0;
}

LRESULT MainWindow::HandleMessage(UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	switch (uMsg)
	{
	case WM_CREATE:
		if (FAILED(D2D1CreateFactory(
			D2D1_FACTORY_TYPE_SINGLE_THREADED, &pFactory)))
		{
			return -1;  // Fail CreateWindowEx.
		}

		M.load();
		M.init();

		return 0;

	case WM_DESTROY:
		DiscardGraphicsResources();
		SafeRelease(&pFactory);
		PostQuitMessage(0);
		return 0;

	case WM_PAINT:
		OnPaint();
		return 0;

	case WM_SIZE:
		Resize();
		return 0;
	}
	return DefWindowProc(m_hwnd, uMsg, wParam, lParam);
}


