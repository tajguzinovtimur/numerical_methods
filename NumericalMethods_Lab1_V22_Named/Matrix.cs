using System.Text.Json;

public static class Matrix
{
    public static double[,] Read(JsonElement element)
    {
        int rows = element.GetArrayLength();
        if (rows == 0) throw new ArgumentException("Пустая матрица");
        int cols = element[0].GetArrayLength();
        if (cols == 0) throw new ArgumentException("Пустая строка матрицы");
        double[,] a = new double[rows, cols];
        for (int i = 0; i < rows; i++)
        {
            if (element[i].GetArrayLength() != cols) throw new ArgumentException("Разная длина строк матрицы");
            for (int j = 0; j < cols; j++) a[i, j] = element[i][j].GetDouble();
        }
        return a;
    }

    public static double[] Vector(JsonElement element)
    {
        double[] a = new double[element.GetArrayLength()];
        for (int i = 0; i < a.Length; i++) a[i] = element[i].GetDouble();
        return a;
    }

    public static void Square(double[,] a)
    {
        if (a.GetLength(0) == 0 || a.GetLength(0) != a.GetLength(1))
            throw new ArgumentException("Матрица должна быть квадратной");
    }

    public static double[,] Identity(int n)
    {
        double[,] a = new double[n, n];
        for (int i = 0; i < n; i++) a[i, i] = 1;
        return a;
    }

    public static double[,] Multiply(double[,] a, double[,] b)
    {
        int rows = a.GetLength(0), inner = a.GetLength(1), cols = b.GetLength(1);
        if (inner != b.GetLength(0)) throw new ArgumentException("Несовместимые матрицы");
        double[,] result = new double[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int k = 0; k < inner; k++)
                for (int j = 0; j < cols; j++) result[i, j] += a[i, k] * b[k, j];
        return result;
    }

    public static double[] Multiply(double[,] a, double[] x)
    {
        if (a.GetLength(1) != x.Length) throw new ArgumentException("Несовместимые размеры");
        double[] result = new double[a.GetLength(0)];
        for (int i = 0; i < result.Length; i++)
            for (int j = 0; j < x.Length; j++) result[i] += a[i, j] * x[j];
        return result;
    }

    public static double[,] Transpose(double[,] a)
    {
        double[,] t = new double[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) t[j, i] = a[i, j];
        return t;
    }

    public static double Difference(double[,] a, double[,] b)
    {
        if (a.GetLength(0) != b.GetLength(0) || a.GetLength(1) != b.GetLength(1))
            throw new ArgumentException("Несовместимые размеры");
        double result = 0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                result = Math.Max(result, Math.Abs(a[i, j] - b[i, j]));
        return result;
    }

    public static double Difference(double[] a, double[] b)
    {
        double result = 0;
        for (int i = 0; i < a.Length; i++) result = Math.Max(result, Math.Abs(a[i] - b[i]));
        return result;
    }

    public static void Print(string name, double[,] a)
    {
        Console.WriteLine($"\n{name}:");
        for (int i = 0; i < a.GetLength(0); i++)
        {
            for (int j = 0; j < a.GetLength(1); j++) Console.Write($"{a[i, j],15:F9}");
            Console.WriteLine();
        }
    }

    public static void Print(string name, double[] a)
    {
        Console.WriteLine($"\n{name}:");
        for (int i = 0; i < a.Length; i++) Console.WriteLine($"{i + 1}: {a[i]:F10}");
    }
}
