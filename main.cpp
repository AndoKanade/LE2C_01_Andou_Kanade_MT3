#include <Novice.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include <assert.h>
const char kWindowTitle[] = "LE2C_01_アンドウ_カナデ_MT3";

typedef struct Matrix4x4 {

  float m[4][4];
} Matrix4x4;

Matrix4x4 Add(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result;
  result.m[0][0] = m1.m[0][0] + m2.m[0][0];
  result.m[0][1] = m1.m[0][1] + m2.m[0][1];
  result.m[0][2] = m1.m[0][2] + m2.m[0][2];
  result.m[0][3] = m1.m[0][3] + m2.m[0][3];
  result.m[1][0] = m1.m[1][0] + m2.m[1][0];
  result.m[1][1] = m1.m[1][1] + m2.m[1][1];
  result.m[1][2] = m1.m[1][2] + m2.m[1][2];
  result.m[1][3] = m1.m[1][3] + m2.m[1][3];
  result.m[2][0] = m1.m[2][0] + m2.m[2][0];
  result.m[2][1] = m1.m[2][1] + m2.m[2][1];
  result.m[2][2] = m1.m[2][2] + m2.m[2][2];
  result.m[2][3] = m1.m[2][3] + m2.m[2][3];
  result.m[3][0] = m1.m[3][0] + m2.m[3][0];
  result.m[3][1] = m1.m[3][1] + m2.m[3][1];
  result.m[3][2] = m1.m[3][2] + m2.m[3][2];
  result.m[3][3] = m1.m[3][3] + m2.m[3][3];
  return result;
}

Matrix4x4 Subtract(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result;
  result.m[0][0] = m1.m[0][0] - m2.m[0][0];
  result.m[0][1] = m1.m[0][1] - m2.m[0][1];
  result.m[0][2] = m1.m[0][2] - m2.m[0][2];
  result.m[0][3] = m1.m[0][3] - m2.m[0][3];
  result.m[1][0] = m1.m[1][0] - m2.m[1][0];
  result.m[1][1] = m1.m[1][1] - m2.m[1][1];
  result.m[1][2] = m1.m[1][2] - m2.m[1][2];
  result.m[1][3] = m1.m[1][3] - m2.m[1][3];
  result.m[2][0] = m1.m[2][0] - m2.m[2][0];
  result.m[2][1] = m1.m[2][1] - m2.m[2][1];
  result.m[2][2] = m1.m[2][2] - m2.m[2][2];
  result.m[2][3] = m1.m[2][3] - m2.m[2][3];
  result.m[3][0] = m1.m[3][0] - m2.m[3][0];
  result.m[3][1] = m1.m[3][1] - m2.m[3][1];
  result.m[3][2] = m1.m[3][2] - m2.m[3][2];
  result.m[3][3] = m1.m[3][3] - m2.m[3][3];
  return result;
}

Matrix4x4 Multiply(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result;
  result.m[0][0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] +
                   m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0];
  result.m[0][1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] +
                   m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1];
  result.m[0][2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] +
                   m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2];
  result.m[0][3] = m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] +
                   m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3];
  result.m[1][0] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] +
                   m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0];
  result.m[1][1] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] +
                   m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1];
  result.m[1][2] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] +
                   m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2];
  result.m[1][3] = m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] +
                   m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3];
  result.m[2][0] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] +
                   m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0];
  result.m[2][1] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] +
                   m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1];
  result.m[2][2] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] +
                   m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2];
  result.m[2][3] = m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] +
                   m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3];
  result.m[3][0] = m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] +
                   m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0];
  result.m[3][1] = m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] +
                   m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1];
  result.m[3][2] = m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] +
                   m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2];
  result.m[3][3] = m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] +
                   m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3];
  return result;
}

Matrix4x4 Inverse(const Matrix4x4 &m) {
  Matrix4x4 result;

  // 行列の行列式を計算
  float det =
      m.m[0][0] *
          (m.m[1][1] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) -
           m.m[1][2] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) +
           m.m[1][3] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1])) -
      m.m[0][1] *
          (m.m[1][0] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) -
           m.m[1][2] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) +
           m.m[1][3] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0])) +
      m.m[0][2] *
          (m.m[1][0] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) -
           m.m[1][1] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) +
           m.m[1][3] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) -
      m.m[0][3] * (m.m[1][0] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1]) -
                   m.m[1][1] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0]) +
                   m.m[1][2] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0]));

  if (det == 0) {
    // 行列式がゼロの場合、逆行列は存在しません
    return result; // 逆行列は存在しないのでゼロ行列を返す
  }

  // 行列式の逆数を計算
  float invDet = 1.0f / det;

  // 各要素を余因子行列から計算
  result.m[0][0] =
      (m.m[1][1] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) -
       m.m[1][2] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) +
       m.m[1][3] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1])) *
      invDet;
  result.m[0][1] =
      (-m.m[0][1] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) +
       m.m[0][2] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) -
       m.m[0][3] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1])) *
      invDet;
  result.m[0][2] =
      (m.m[0][1] * (m.m[1][2] * m.m[3][3] - m.m[1][3] * m.m[3][2]) -
       m.m[0][2] * (m.m[1][1] * m.m[3][3] - m.m[1][3] * m.m[3][1]) +
       m.m[0][3] * (m.m[1][1] * m.m[3][2] - m.m[1][2] * m.m[3][1])) *
      invDet;
  result.m[0][3] =
      (-m.m[0][1] * (m.m[1][2] * m.m[2][3] - m.m[1][3] * m.m[2][2]) +
       m.m[0][2] * (m.m[1][1] * m.m[2][3] - m.m[1][3] * m.m[2][1]) -
       m.m[0][3] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1])) *
      invDet;

  result.m[1][0] =
      (-m.m[1][0] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) +
       m.m[1][2] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) -
       m.m[1][3] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0])) *
      invDet;
  result.m[1][1] =
      (m.m[0][0] * (m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2]) -
       m.m[0][2] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) +
       m.m[0][3] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0])) *
      invDet;
  result.m[1][2] =
      -(m.m[0][0] * (m.m[1][2] * m.m[3][3] - m.m[1][3] * m.m[3][2]) -
       m.m[0][2] * (m.m[1][0] * m.m[3][3] - m.m[1][3] * m.m[3][0]) +
       m.m[0][3] * (m.m[1][0] * m.m[3][2] - m.m[1][2] * m.m[3][0])) *
      invDet;

  result.m[1][3] =
      (m.m[0][0] * (m.m[1][2] * m.m[2][3] - m.m[1][3] * m.m[2][2]) -
       m.m[0][2] * (m.m[1][0] * m.m[2][3] - m.m[1][3] * m.m[2][0]) +
       m.m[0][3] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0])) *
      invDet;

    result.m[2][0] =
      (m.m[1][0] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) -
       m.m[1][1] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) +
       m.m[1][3] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) *
      invDet;

  result.m[2][1] =
      (-m.m[0][0] * (m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1]) +
       m.m[0][1] * (m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0]) -
       m.m[0][3] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) *
      invDet;

  result.m[2][2] =
      (m.m[0][0] * (m.m[1][1] * m.m[3][3] - m.m[1][3] * m.m[3][1]) -
       m.m[0][1] * (m.m[1][0] * m.m[3][3] - m.m[1][3] * m.m[3][0]) +
       m.m[0][3] * (m.m[1][0] * m.m[3][1] - m.m[1][1] * m.m[3][0])) *
      invDet;

  result.m[2][3] =
      (-m.m[0][0] * (m.m[1][1] * m.m[2][3] - m.m[1][3] * m.m[2][1]) +
       m.m[0][1] * (m.m[1][0] * m.m[2][3] - m.m[1][3] * m.m[2][0]) -
       m.m[0][3] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0])) *
      invDet;

  result.m[3][0] =
      (-m.m[1][0] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1]) +
       m.m[1][1] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0]) -
       m.m[1][2] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) *
      invDet;

  result.m[3][1] =
      (m.m[0][0] * (m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1]) -
       m.m[0][1] * (m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0]) +
       m.m[0][2] * (m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0])) *
      invDet;

  result.m[3][2] =
      (-m.m[0][0] * (m.m[1][1] * m.m[3][2] - m.m[1][2] * m.m[3][1]) +
       m.m[0][1] * (m.m[1][0] * m.m[3][2] - m.m[1][2] * m.m[3][0]) -
       m.m[0][2] * (m.m[1][0] * m.m[3][1] - m.m[1][1] * m.m[3][0])) *
      invDet;

result.m[3][3] =
      (m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
       m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
       m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0])) *
      invDet;

  return result;
}

Matrix4x4 Transpose(const Matrix4x4 &m) {
  Matrix4x4 result;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      result.m[i][j] = m.m[j][i];
    }
  }
  return result;
}

Matrix4x4 MakeIdentity4x4() {
  Matrix4x4 result;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j) {
        result.m[i][j] = 1.0f;
      } else {
        result.m[i][j] = 0.0f;
      }
    }
  }
  return result;
}

static const int kRowheight = 20;
static const int kColumnWidth = 60;

void Matrix4x4ScreenPrintf(int x, int y, const Matrix4x4 &matrix,
                           const char *lavel) {

  Novice::ScreenPrintf(x, y, lavel);
  for (int row = 0; row < 4; row++) {
    for (int column = 0; column < 4; column++) {
      Novice::ScreenPrintf(x + column * kColumnWidth,
                           y + (row + 1) * kRowheight, "%6.02f",
                           matrix.m[row][column]);
    }
  }
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, 1280, 720);

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

  Matrix4x4 m1 = {3.2f, 0.7f, 9.6f, 4.4f, 5.5f, 1.3f, 7.8f, 2.1f,
                  6.9f, 8.0f, 2.6f, 1.0f, 0.5f, 7.2f, 5.1f, 3.3f};

  Matrix4x4 m2 = {4.1f, 6.5f, 3.3f, 2.2f, 8.8f, 0.6f, 9.9f, 7.7f,
                  1.1f, 5.5f, 6.6f, 0.0f, 3.3f, 9.9f, 8.8f, 2.2f};

  // ウィンドウの×ボタンが押されるまでループ
  while (Novice::ProcessMessage() == 0) {
    // フレームの開始
    Novice::BeginFrame();

    // キー入力を受け取る
    memcpy(preKeys, keys, 256);
    Novice::GetHitKeyStateAll(keys);

    ///
    /// ↓更新処理ここから
    ///

    Matrix4x4 resultAdd = Add(m1, m2);
    Matrix4x4 resultSubtract = Subtract(m1, m2);
    Matrix4x4 resultMultiply = Multiply(m1, m2);
    Matrix4x4 resultInverseM1 = Inverse(m1);
    Matrix4x4 resultInverseM2 = Inverse(m2);
    Matrix4x4 transposeM1 = Transpose(m1);
    Matrix4x4 transposeM2 = Transpose(m2);
    Matrix4x4 identity = MakeIdentity4x4();

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    Matrix4x4ScreenPrintf(0, 0, resultAdd, "Add");
    Matrix4x4ScreenPrintf(0, kRowheight * 5, resultSubtract, "Subtract");
    Matrix4x4ScreenPrintf(0, kRowheight * 10, resultMultiply, "Multiply");
    Matrix4x4ScreenPrintf(0, kRowheight * 15, resultInverseM1, "InverseM1");
    Matrix4x4ScreenPrintf(0, kRowheight * 20, resultInverseM2, "InverseM2");
    Matrix4x4ScreenPrintf(kColumnWidth * 5, 0, transposeM1, "TransposeM1");
    Matrix4x4ScreenPrintf(kColumnWidth * 5, kRowheight * 5, transposeM2,
                          "TransposeM2");
    Matrix4x4ScreenPrintf(kColumnWidth * 5, kRowheight * 10, identity,
                          "Identity");

    ///
    /// ↑描画処理ここまで
    ///

    // フレームの終了
    Novice::EndFrame();

    // ESCキーが押されたらループを抜ける
    if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
      break;
    }
  }

  // ライブラリの終了
  Novice::Finalize();
  return 0;
}
