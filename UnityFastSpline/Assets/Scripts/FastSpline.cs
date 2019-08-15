using System.Collections;
using System.Collections.Generic;
using UnityEngine;

#if UNITY_EDITOR
using UnityEditor;

[CustomEditor(typeof(FastSpline))]
public class FastSplineEditor
    : Editor
{
    private int selectedIndex = -1;

    private bool preview;
    private float previewTime;

    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();

        var self = target as FastSpline;

        preview = GUILayout.Toggle(preview, "Preview");
        previewTime = EditorGUILayout.Slider("Preview Time", previewTime, 0.0f, 1.0f);
    }

    public void OnSceneGUI()
    {
        var self = target as FastSpline;
        var steps = 10;

        for (int i = 1; i < steps; ++i)
        {
            var t0 = 1.0f / (float)steps * (float)(i - 1);
            var t1 = 1.0f / (float)steps * (float)(i);
            var p0 = self.CalculatePosition(t0);
            var p1 = self.CalculatePosition(t1);
            var w0 = self.transform.TransformPoint(p0);
            var w1 = self.transform.TransformPoint(p1);

            Handles.DrawLine(w0, w1);
        }

        for (int i = 0; i < self.points.Count; ++i)
        {
            var p0 = self.points[i];
            var w0 = self.transform.TransformPoint(p0);

            if (selectedIndex == i)
            {
                Handles.color = Color.blue;
                Handles.Button(w0, Quaternion.identity, 0.5f, 0.5f, Handles.RectangleHandleCap);

                var refPosition = p0;
                var refRotation = Quaternion.identity;

                refPosition = self.transform.TransformPoint(refPosition);

                Handles.TransformHandle(ref refPosition, ref refRotation);

                refPosition = self.transform.InverseTransformPoint(refPosition);

                self.points[i] = refPosition;
            }
            else
            {
                Handles.color = Color.white;

                if (Handles.Button(w0, Quaternion.identity, 0.5f, 0.5f, Handles.RectangleHandleCap))
                    selectedIndex = i;
            }
        }

        if (preview)
        {
            var p0 = self.CalculatePosition(previewTime);
            var w0 = self.transform.TransformPoint(p0);

            Handles.color = Color.red;
            Handles.DrawSphere(-1, w0, Quaternion.identity, 0.5f);
        }
    }
}
#endif


public static class FastSplineMath
{
    public static Vector3 QuadraticInterpolate(List<Vector3> points, int index, float t)
    {
        var count = points.Count;

        var index0 = Mathf.Clamp(index + 0, 0, count - 1);
        var index1 = Mathf.Clamp(index + 1, 0, count - 1);
        var index2 = Mathf.Clamp(index + 2, 0, count - 1);

        var p0 = points[index0];
        var p1 = points[index1];
        var p2 = points[index2];

        return QuadraticInterpolate(p0, p1, p2, t);
    }

    public static Vector3 QuadraticInterpolate(Vector3 p0, Vector3 p1, Vector3 p2, float t)
    {
        var tt = t * t;
        var u = 1.0f - t;
        var uu = u * u;
        var p = p0 * uu + p1 * 2.0f * u * t + p2 * tt;

        return p;
    }

    public static Vector3 CubicInterpolate(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, float t)
    {
        var u = 1.0f - t;
        var tt = t * t;
        var uu = u * u;
        var uuu = uu * u;
        var ttt = tt * t;
        var p = uuu * p0;

        p += 3.0f * uu * t * p1;
        p += 3.0f* u * tt * p2;
        p += ttt * p3;

        return p;
    }
}

public class FastSpline
    : MonoBehaviour
{
    public List<Vector3> points;

    public FastSpline()
    {
        if (points == null)
        {
            points = new List<Vector3>();
            points.Add(Vector3.zero);
            points.Add(Vector3.zero);
            points.Add(Vector3.zero);
        }
    }

    public Vector3 CalculatePosition(float t)
    {
        var index = (int)(t * (float)points.Count);
        var time = (t / (float)points.Count) * 3.0f;
        var point = FastSplineMath.QuadraticInterpolate(points, index, t);

        return point;
    }

#if UNITY_EDITOR
    public void OnDrawGizmos()
    {
        var steps = 10;

        for (int i = 1; i < steps; ++i)
        {
            var t0 = 1.0f / (float)steps * (float)(i - 1);
            var t1 = 1.0f / (float)steps * (float)(i);
            var p0 = CalculatePosition(t0);
            var p1 = CalculatePosition(t1);
            var w0 = transform.TransformPoint(p0);
            var w1 = transform.TransformPoint(p1);

            Gizmos.color = Color.Lerp(Color.blue, Color.clear, 0.5f);
            Gizmos.DrawLine(w0, w1);
        }
    }
#endif
}
