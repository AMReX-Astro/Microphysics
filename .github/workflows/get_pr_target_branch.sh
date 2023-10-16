#!/bin/bash

# Get the PR number from the GitHub context
PR_NUMBER=${{ github.event.pull_request.number }}

# Call the GitHub API to get the PR details
PR_DETAILS=$(curl \
                   -H "Authorization: token ${{ secrets.GITHUB_TOKEN }}" \
                     -H "Accept: application/vnd.github.v3+json" \
                       "https://api.github.com/repos/${{ github.repository }}/pulls/$PR_NUMBER")

# Extract the base ref (target branch) from the PR details
TARGET_BRANCH=$(echo "$PR_DETAILS" | jq -r '.base.ref')

# Output the target branch
echo "::set-output name=target_branch::$TARGET_BRANCH"
